classdef (Abstract) mme_branch_opf_ac < mp.mme_branch_opf
% mp.mme_branch_opf_ac - Math model element abstract base class for branch for AC OPF.
%
% Math model element abstract base class for branch elements, including
% transmission lines and transformers, for AC OPF problems.
%
% Implements methods for adding of branch flow constraints and for updating
% the output data in the corresponding data model element for in-service
% branches from the math model solution.

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

            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);
            ibr = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(ibr);      %% number of constrained branches
            mm.userdata.flow_constrained_branch_idx = ibr;

            if nl2
                %% port indexes
                nl = nme.nk;
                idx = [ibr; nl+ibr];

                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                lim_type = upper(mpopt.opf.flow_lim(1));
                if lim_type == 'S'
                    fcn_flow = @(x)port_apparent_power_lim_fcn(nme, ...
                        mm.convert_x_m2n(x, nm), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_apparent_power_lim_hess(nme, ...
                        mm.convert_x_m2n(x, nm), lam, nm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)port_active_power_lim_fcn(nme, ...
                        mm.convert_x_m2n(x, nm), nm, idx, [flow_max; flow_max]);
                    hess_flow = @(x, lam)port_active_power_lim_hess(nme, ...
                        mm.convert_x_m2n(x, nm), lam, nm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)port_active_power2_lim_fcn(nme, ...
                        mm.convert_x_m2n(x, nm), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_active_power2_lim_hess(nme, ...
                        mm.convert_x_m2n(x, nm), lam, nm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)port_current_lim_fcn(nme, ...
                        mm.convert_x_m2n(x, nm), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_current_lim_hess(nme, ...
                        mm.convert_x_m2n(x, nm), lam, nm, idx);
                else
                    error('mp.mme_branch_opf_ac.add_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end

                mm.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% branch active power flow
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1));
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2));

            %% shadow prices on branch flow constraints
            ibr = mm.userdata.flow_constrained_branch_idx;
            mu_flow_fr_ub = zeros(nme.nk, 1);
            mu_flow_to_ub = mu_flow_fr_ub;
            if length(ibr)
                lim_type = upper(mpopt.opf.flow_lim(1));
                nni = mm.get_idx('nli');
                lambda = mm.soln.lambda;
                if lim_type == 'P'
                    mu_flow_fr_ub(ibr) = lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf);
                    mu_flow_to_ub(ibr) = lambda.ineqnonlin(nni.i1.St:nni.iN.St);
                else
                    rate_a = dme.rate_a(ibr);
                    mu_flow_fr_ub(ibr) = 2 * lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf) .* rate_a;
                    mu_flow_to_ub(ibr) = 2 * lambda.ineqnonlin(nni.i1.St:nni.iN.St) .* rate_a;
                end
            end

            %% shadow prices on angle difference limits
            [mu_vad_lb, mu_vad_ub] = obj.ang_diff_prices(mm, nme);

            %% update in the data model
            dme.tab.pl_fr(dme.on) = real(S_fr) * dm.base_mva;
            dme.tab.ql_fr(dme.on) = imag(S_fr) * dm.base_mva;
            dme.tab.pl_to(dme.on) = real(S_to) * dm.base_mva;
            dme.tab.ql_to(dme.on) = imag(S_to) * dm.base_mva;
            dme.tab.mu_flow_fr_ub(dme.on) = mu_flow_fr_ub / dm.base_mva;
            dme.tab.mu_flow_to_ub(dme.on) = mu_flow_to_ub / dm.base_mva;
            dme.tab.mu_vad_lb(dme.on) = mu_vad_lb * pi/180;
            dme.tab.mu_vad_ub(dme.on) = mu_vad_ub * pi/180;
        end
    end     %% methods
end         %% classdef
