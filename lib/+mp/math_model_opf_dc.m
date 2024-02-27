classdef math_model_opf_dc < mp.math_model_opf
% mp.math_model_opf_dc - Optimal Power flow (OPF) **math model** for DC formulation.
%
% Provides formulation-specific and OPF-specific subclasses for elements.
%
% Provides implementation of nodal balance constraint method and setup
% of solver options.
%
% Implements convert_x_m2n() to convert from math model state to network
% model state.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = math_model_opf_dc()
            %

            obj@mp.math_model_opf();
            obj.element_classes = { @mp.mme_bus_opf_dc, @mp.mme_gen_opf_dc, ...
                @mp.mme_load_pf_dc, @mp.mme_branch_opf_dc, @mp.mme_shunt_pf_dc };
        end

        function tag = form_tag(obj)
            %

            tag = 'dc';
        end

        function name = form_name(obj)
            %

            name = 'DC';
        end

        function [vx, z, x] = convert_x_m2n(obj, mmx, nm)
            %

            nm_vars = obj.update_nm_vars(mmx, nm);

            %% convert (real) math model x to network model x
            vx = nm_vars.va;
            z  = nm_vars.z;
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            [B, K, p] = nm.get_params();

            %% power balance constraints
            C = nm.C;
            Amis = C * [B*C' K*nm.D'];
            bmis = -C * p;
            obj.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                [nm.va.order nm.z.order]);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %

            opt = mpopt2qpopt(mpopt, obj.problem_type());

            switch opt.alg
                case {'MIPS', 'IPOPT'}
                    if mpopt.opf.start < 2      %% initialize interior point
                        opt.x0 = obj.interior_x0(obj, nm, dm);
                    end
                case 'OSQP'
                    opt.x0 = [];        %% disable provided starting point
            end
        end
    end     %% methods
end         %% classdef
