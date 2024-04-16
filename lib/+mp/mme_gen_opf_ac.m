classdef mme_gen_opf_ac < mp.mme_gen_opf
% mp.mme_gen_opf_ac - Math model element for generator for AC OPF.
%
% Math model element class for generator elements for AC OPF problems.
%
% Implements methods for buliding and adding PQ capability constraints,
% dispatchable load power factor constraints, polynomial costs, and for
% updating the output data in the corresponding data model element for
% in-service generators from the math model solution.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %

            %% generator PQ capability curve constraints
            [Apqh, ubpqh, Apql, ubpql, Apqdata] = ...
                obj.pq_capability_constraint(dm.elements.gen, dm.base_mva);
            mm.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
            mm.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
            mm.userdata.Apqdata = Apqdata;

            %% dispatchable load constant power factor constraint
            [Avl, lvl, uvl] = obj.disp_load_constant_pf_constraint(dm);
            if ~isempty(Avl)
                mm.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});    %% nvl
            end

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'Qg', 'y'});
            end

            %% call parent
            add_constraints@mp.mme_gen_opf(obj, mm, nm, dm, mpopt);
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            %

            %% call parent
            add_costs@mp.mme_gen_opf(obj, mm, nm, dm, mpopt);

            %% costs on reactive dispatch
            if ~isempty(obj.cost.poly_q)
                %% (quadratic) polynomial costs on Qg
                if obj.cost.poly_q.have_quad_cost
                    mm.add_quad_cost('polQg', obj.cost.poly_q.Q, obj.cost.poly_q.c, obj.cost.poly_q.k, {'Qg'});
                end

                %% (order 3 and higher) polynomial costs on Qg
                if ~isempty(obj.cost.poly_q.i3)
                    dme = obj.data_model_element(dm);
                    cost_Qg = @(xx)mp.cost_table.poly_cost_fcn( ...
                        xx, dm.base_mva, ...
                        dme.tab.cost_qg.poly_coef(dme.on, :), ...
                        obj.cost.poly_q.i3);
                    mm.add_nln_cost('polQg', 1, cost_Qg, {'Qg'});
                end
            end
        end

        function [Ah, uh, Al, ul, data] = pq_capability_constraint(obj, dme, base_mva)
            % from legacy :func:`makeApq`

            gen = dme.tab(dme.on, :);

            %% data dimensions
            ng = size(gen, 1);      %% number of dispatchable injections

            %% which generators require additional linear constraints
            %% (in addition to simple box constraints) on (pg,qg) to correctly
            %% model their PQ capability curves
            ipqh = find( obj.has_pq_cap(gen, 'U') );
            ipql = find( obj.has_pq_cap(gen, 'L') );
            npqh = size(ipqh, 1);   %% number of general PQ capability curves (upper)
            npql = size(ipql, 1);   %% number of general PQ capability curves (lower)

            %% make Ah if there is a need to add general PQ capability curves;
            %% use normalized coefficient rows so multipliers have right scaling
            %% in $$/pu
            if npqh > 0
                data.h = [gen.qc1_ub(ipqh)-gen.qc2_ub(ipqh), gen.pc2(ipqh)-gen.pc1(ipqh)];
                uh = data.h(:, 1) .* gen.pc1(ipqh) + data.h(:, 2) .* gen.qc1_ub(ipqh);
                for i=1:npqh,
                    tmp = norm(data.h(i,:));
                    data.h(i,:) = data.h(i, :) / tmp;
                    uh(i) = uh(i) / tmp;
                end
                Ah = sparse([1:npqh, 1:npqh]', [ipqh; ipqh+ng], ...
                            data.h(:), npqh, 2*ng);
                uh = uh / base_mva;
            else
                data.h = [];
                Ah  = sparse(0, 2*ng);
                uh = [];
            end

            %% similarly Al
            if npql > 0
                data.l = [gen.qc2_lb(ipql)-gen.qc1_lb(ipql), gen.pc1(ipql)-gen.pc2(ipql)];
                ul= data.l(:, 1) .* gen.pc1(ipql) + data.l(:, 2) .* gen.qc1_lb(ipql) ;
                for i=1:npql,
                    tmp = norm(data.l(i, : ));
                    data.l(i, :) = data.l(i, :) / tmp;
                    ul(i) = ul(i) / tmp;
                end
                Al = sparse([1:npql, 1:npql]', [ipql; ipql+ng], ...
                            data.l(:), npql, 2*ng);
                ul = ul / base_mva;
            else
                data.l = [];
                Al  = sparse(0, 2*ng);
                ul = [];
            end

            data.ipql = ipql;
            data.ipqh = ipqh;
        end

        function TorF = has_pq_cap(obj, gen, upper_lower)
            % from legacy :func:`hasPQcap`

            %% default value
            if nargin < 3
                upper_lower = 'B';  %% look at both top and bottom by default
            end

            %% for which gens is it specified
            k = find( gen.pc1 | gen.pc2 );
            ng = size(gen, 1);

            if isempty(k)
                TorF = zeros(ng, 1);
            else
                %% eliminate cases where QMIN = QMAX = QC
                kk = find(  gen.qg_lb(k) == gen.qg_ub(k) & ...
                            gen.qg_lb(k) == gen.qc1_ub(k) & ...
                            gen.qg_lb(k) == gen.qc1_lb(k) & ...
                            gen.qg_lb(k) == gen.qc2_ub(k) & ...
                            gen.qg_lb(k) == gen.qc2_lb(k) );
                k(kk) = [];

                %% check for errors in capability curve data
                if any( gen.pc1(k) >= gen.pc2(k) )
                    error('mp.mme_gen_opf_ac.has_pq_cap: must have pc1 < pc2');
                end
                if any( gen.qc2_ub(k) <= gen.qc2_lb(k) & gen.qc1_ub(k) <= gen.qc1_lb(k) )
                    error('mp.mme_gen_opf_ac.has_pq_cap: capability curve defines an empty set');
                end

                %% for which gens is it specified
                k = find( gen.pc1 ~= gen.pc2 );
                L = zeros(ng, 1);
                U = zeros(ng, 1);
                dPc = gen.pc2(k) - gen.pc1(k);

                if ~strcmp(upper_lower, 'U')    %% include lower constraint
                    dQc = gen.qc2_lb(k) - gen.qc1_lb(k);
                    Qmin_at_Pmin = gen.qc1_lb(k) + (gen.pg_lb(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    Qmin_at_Pmax = gen.qc1_lb(k) + (gen.pg_ub(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    L(k) = Qmin_at_Pmin > gen.qg_lb(k) | Qmin_at_Pmax > gen.qg_lb(k);
                end

                if ~strcmp(upper_lower, 'L')    %% include upper constraint
                    dQc = gen.qc2_ub(k) - gen.qc1_ub(k);
                    Qmax_at_Pmin = gen.qc1_ub(k) + (gen.pg_lb(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    Qmax_at_Pmax = gen.qc1_ub(k) + (gen.pg_ub(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    U(k) = Qmax_at_Pmin < gen.qg_ub(k) | Qmax_at_Pmax < gen.qg_ub(k);
                end

                TorF = L | U;
            end
        end

        function [A, l, u] = disp_load_constant_pf_constraint(obj, dm)
            % from legacy :func:`makeAvl`

            dme = obj.data_model_element(dm);

            %% data dimensions
            ng = dme.n;     %% number of dispatchable injections
            pg = dme.pg_start;
            qg = dme.qg_start;
            pg_lb = dme.pg_lb;
            qg_lb = dme.qg_lb;
            qg_ub = dme.qg_ub;

            ivl = find( dme.isload(dme.on) & (qg_lb ~= 0 | qg_ub ~= 0) );
            nvl  = size(ivl, 1);  %% number of dispatchable loads

            %% at least one of the Q limits must be zero (corresponding to pg_ub == 0)
            if any( qg_lb(ivl) ~= 0 & qg_ub(ivl) ~= 0 )
                k = find(qg_lb(ivl) ~= 0 & qg_ub(ivl) ~= 0);
                gidx = dme.tab.uid(dme.on(ivl(k)));
                s = sprintf('Invalid qg limits for dispatchable load in row %d of gen table\n', gidx);
                error('mp.mme_gen_opf_ac.disp_load_constant_pf_constraint: Either qg_lb or qg_ub must be equal to zero for each dispatchable load.\n%s', s);
            end

            %% Initial values of PG and QG must be consistent with specified
            %% power factor
            Qlim = (qg_lb(ivl) == 0) .* qg_ub(ivl) + ...
                (qg_ub(ivl) == 0) .* qg_lb(ivl);
            if any( abs( qg(ivl) - pg(ivl) .* Qlim ./ pg_lb(ivl) ) > 1e-6 )
                k = find(abs( qg(ivl) - pg(ivl) .* Qlim ./ pg_lb(ivl) ) > 1e-6);
                gidx = dme.tab.uid(dme.on(ivl(k)));
                s = sprintf('qg for dispatchable load in row %d of gen table must be pg * %g\n', [gidx Qlim ./ pg_lb(ivl)]');
                error('mp.mme_gen_opf_ac.disp_load_constant_pf_constraint: %s\n         %s\n         %s\n         %s\n%s', ...
                    'For a dispatchable load, pg and qg must be consistent', ...
                    'with the power factor defined by pg_lb and the relevant', ...
                    '(non-zero) qg_lb or qg_ub limit.', ...
                    'Note: Setting pg = qg = 0 satisfies this condition.', s);
            end

            %% make A, l, u, for l <= A * [pg; qg] <= u
            if nvl > 0
                xx = pg_lb(ivl);
                yy = Qlim;
                pftheta = atan2(yy, xx);
                pc = sin(pftheta);
                qc = -cos(pftheta);
                ii = [ (1:nvl)'; (1:nvl)' ];
                jj = [ ivl; ivl+ng ];
                A = sparse(ii, jj, [pc; qc], nvl, 2*ng);
                l = zeros(nvl, 1);
                u = l;
            else
                A = sparse(0, 2*ng);
                l = [];
                u = [];
            end
        end

        function build_cost_params(obj, dm)
            %
            dme = obj.data_model_element(dm);
            obj.cost = dme.build_cost_params(dm, 0);
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% generator active power
            ss = nm.get_idx('state');
            Sg = nm.soln.z(ss.i1.gen:ss.iN.gen);
            vm_setpoint = abs(nme.C' * nm.soln.v);

            %% shadow prices on generator limits
            [vv, ll] = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_pg_lb = lambda.lower(vv.i1.Pg:vv.iN.Pg);
            mu_pg_ub = lambda.upper(vv.i1.Pg:vv.iN.Pg);
            mu_qg_lb = lambda.lower(vv.i1.Qg:vv.iN.Qg);
            mu_qg_ub = lambda.upper(vv.i1.Qg:vv.iN.Qg);

            %% gen PQ capability curve multipliers - based on update_mupq()
            if ll.N.PQh > 0 || ll.N.PQl > 0
                d = mm.get_userdata('Apqdata');

                %% combine original limit multipliers into single value
                muP = mu_pg_ub - mu_pg_lb;
                muQ = mu_qg_ub - mu_qg_lb;

                %% add P and Q components of multipliers on upper sloped constraint
                if ~isempty(d.ipqh)
                    mu_PQh = lambda.mu_l(ll.i1.PQh:ll.iN.PQh) - lambda.mu_u(ll.i1.PQh:ll.iN.PQh);
                    muP(d.ipqh) = muP(d.ipqh) - mu_PQh .* d.h(:,1);
                    muQ(d.ipqh) = muQ(d.ipqh) - mu_PQh .* d.h(:,2);
                end

                %% add P and Q components of multipliers on lower sloped constraint
                if ~isempty(d.ipql)
                    mu_PQl = lambda.mu_l(ll.i1.PQl:ll.iN.PQl) - lambda.mu_u(ll.i1.PQl:ll.iN.PQl);
                    muP(d.ipql) = muP(d.ipql) - mu_PQl .* d.l(:,1);
                    muQ(d.ipql) = muQ(d.ipql) - mu_PQl .* d.l(:,2);
                end

                %% split back into upper and lower multipliers based on sign
                mu_pg_ub = (muP > 0) .*  muP;
                mu_pg_lb = (muP < 0) .* -muP;
                mu_qg_ub = (muQ > 0) .*  muQ;
                mu_qg_lb = (muQ < 0) .* -muQ;
            end

            %% update in the data model
            dme.tab.pg(dme.on) = real(Sg) * dm.base_mva;
            dme.tab.qg(dme.on) = imag(Sg) * dm.base_mva;
            dme.tab.vm_setpoint(dme.on) = vm_setpoint;
            dme.tab.mu_pg_lb(dme.on) = mu_pg_lb / dm.base_mva;
            dme.tab.mu_pg_ub(dme.on) = mu_pg_ub / dm.base_mva;
            dme.tab.mu_qg_lb(dme.on) = mu_qg_lb / dm.base_mva;
            dme.tab.mu_qg_ub(dme.on) = mu_qg_ub / dm.base_mva;
        end
    end     %% methods
end         %% classdef
