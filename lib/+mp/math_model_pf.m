classdef (Abstract) math_model_pf < mp.math_model
% mp.math_model_pf - Abstract base class for power flow (PF) **math model** objects.
%
% Implements setting up of solver options from |MATPOWER| options struct.

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
        function tag = task_tag(obj)
            %

            tag = 'pf';
        end

        function name = task_name(obj)
            %

            name = 'Power Flow';
        end

        function obj = add_costs(obj, nm, dm, mpopt)
            %

            %% no costs for PF and CPF
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
            %

            %% implementation in method in mp.mm_shared_pfcpf class hierarchy
            obj.add_system_vars_pf(nm, dm, mpopt);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %

            %% for AC power flow only, must override for mp.math_model_pf_dc
            switch mpopt.pf.alg
                case 'DEFAULT'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'DEFAULT');
                case {'NR', 'NR-SP', 'NR-SC', 'NR-SH', 'NR-IP', 'NR-IC', 'NR-IH'}
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'NEWTON');
                case {'FDXB', 'FDBX'}
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'FD');
                    opt.fd_opt.jac_approx_fcn = @()obj.fd_jac_approx(nm, dm, mpopt);
                    opt.fd_opt.labels = {'P', 'Q'};
                case 'FSOLVE'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'FSOLVE');
                case 'GS'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'GS');
                    opt.gs_opt.x_update_fcn = ...
                        @(x, f)obj.gs_x_update(x, f, nm, dm, mpopt);
                case 'ZG'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'ZG');
                    zg_x_update = @(x, f)obj.zg_x_update(x, f, nm, dm, mpopt);
                    opt.core_sp = struct(...
                        'alg',              'ZG', ...
                        'name',             'Implicit Z-bus Gauss', ...
                        'default_max_it',   1000, ...
                        'need_jac',         0, ...
                        'update_fcn',       zg_x_update  );
                otherwise
                    error('mp.math_model_pf.solve_opts: invalid value for MPOPT.PF.ALG (%s)', mpopt.pf.alg);
            end
            opt.verbose = mpopt.verbose;
        end
    end     %% methods
end         %% classdef
