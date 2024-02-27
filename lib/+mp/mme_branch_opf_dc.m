classdef mme_branch_opf_dc < mp.mme_branch_opf
% mp.mme_branch_opf_dc - Math model element for branch for DC OPF.
%
% Math model element class for branch elements, including transmission lines
% and transformers, for DC OPF problems.
%
% Implements methods for adding of branch flow and angle difference
% constraints and for updating the output data in the corresponding data
% model element for in-service branches from the math model solution.

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
                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                [B, K, p] = nme.get_params(ibr);
                Af = B * nme.C';
                mm.add_lin_constraint('Pf', Af, -p-flow_max, -p+flow_max, ...
                    nm.va.order);
            end

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = obj.ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            if length(iang)
                mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% branch active power flow
            pp = nm.get_idx('port');
            pl_fr = nm.soln.gp(pp.i1.branch(1):pp.iN.branch(1));
            pl_to = nm.soln.gp(pp.i1.branch(2):pp.iN.branch(2));

            %% shadow prices on branch flow constraints
            ibr = mm.userdata.flow_constrained_branch_idx;
            mu_flow_fr_ub = zeros(nme.nk, 1);
            mu_flow_to_ub = mu_flow_fr_ub;
            if length(ibr)
                ll = mm.get_idx('lin');
                lambda = mm.soln.lambda;
                mu_flow_fr_ub(ibr) = lambda.mu_u(ll.i1.Pf:ll.iN.Pf);
                mu_flow_to_ub(ibr) = lambda.mu_l(ll.i1.Pf:ll.iN.Pf);
            end

            %% shadow prices on angle difference limits
            [mu_vad_lb, mu_vad_ub] = obj.ang_diff_prices(mm, nme);

            %% update in the data model
            dme.tab.pl_fr(dme.on) = pl_fr * dm.base_mva;
            dme.tab.pl_to(dme.on) = pl_to * dm.base_mva;
            dme.tab.mu_flow_fr_ub(dme.on) = mu_flow_fr_ub / dm.base_mva;
            dme.tab.mu_flow_to_ub(dme.on) = mu_flow_to_ub / dm.base_mva;
            dme.tab.mu_vad_lb(dme.on) = mu_vad_lb * pi/180;
            dme.tab.mu_vad_ub(dme.on) = mu_vad_ub * pi/180;
        end
    end     %% methods
end         %% classdef
