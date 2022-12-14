classdef mme_gen_opf_dc < mp.mme_gen_opf

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'y'});
            end

            %% call parent
            add_constraints@mp.mme_gen_opf(obj, mm, nm, dm, mpopt);
        end

        function build_cost_params(obj, dm)
            dme = obj.data_model_element(dm);
            obj.cost = dme.build_cost_params(dm, 1);
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            pg = nm.soln.z(ss.i1.gen:ss.iN.gen) * dm.base_mva;

            %% shadow prices on generator limits
            vv = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_pg_ub = lambda.upper(vv.i1.Pg:vv.iN.Pg);
            mu_pg_lb = lambda.lower(vv.i1.Pg:vv.iN.Pg);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = pg;
            dme.tab.mu_pg_lb(dme.on) = mu_pg_lb / dm.base_mva;
            dme.tab.mu_pg_ub(dme.on) = mu_pg_ub / dm.base_mva;
        end
    end     %% methods
end         %% classdef
