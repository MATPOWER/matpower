classdef math_model_pf_dc < mp.math_model_pf & mp.mm_shared_pfcpf_dc
%MP.MATH_MODEL_PF_DC  MATPOWER mathematical model for DC power flow (PF) problem.
%   ?
%
%   MP.MATH_MODEL_PF_DC ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = math_model_pf_dc()
            obj@mp.math_model_pf();
            obj.element_classes = { @mp.mme_bus_pf_dc, @mp.mme_gen_pf_dc, ...
                @mp.mme_load_pf_dc, @mp.mme_branch_pf_dc, @mp.mme_shunt_pf_dc };
        end

        function tag = form_tag(obj)
            tag = 'dc';
        end

        function name = form_name(obj)
            name = 'DC';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% power balance constraints
            A = ad.B(pvq, pvq);
            b = (ad.Pbus(pvq) - ad.B(pvq, ad.ref) * ad.va(ad.ref));
            obj.add_lin_constraint('Pmis', A, b, b);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %% overrides mp.math_model_pf/solve_opts()
            %% TO DO: move pf.alg to pf.ac.solver and add a
            %%        pf.dc.solver to set the 'leq_opt.solver' option here
            opt = struct( ...
                'verbose',  mpopt.verbose, ...
                'leq_opt',  struct('thresh', 1e5)   );
        end
    end     %% methods
end         %% classdef
