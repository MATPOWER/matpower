classdef dmce_branch_mpc2 < mp.dmc_element % & mp.dmce_branch
% mp.dmce_branch_mpc2 - Data model converter element for branch for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            %
            name = 'branch';
        end

        function df = data_field(obj)
            %
            df = 'branch';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% mapping for each name, default is {'col', []}
            vmap.uid        = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name       = {'cell', ''}; %% empty char
            vmap.status{2}  = BR_STATUS;
            vmap.source_uid = {'cell', ''}; %% empty char
            vmap.bus_fr{2}  = F_BUS;
            vmap.bus_to{2}  = T_BUS;
            vmap.r{2}       = BR_R;
            vmap.x{2}       = BR_X;
            vmap.g_fr       = {'num', 0};   %% zeros
            vmap.b_fr       = {'col', BR_B, 0.5};
            vmap.g_to       = {'num', 0};   %% zeros
            vmap.b_to       = {'col', BR_B, 0.5};
            vmap.sm_ub_a{2} = RATE_A;
            vmap.sm_ub_b{2} = RATE_B;
            vmap.sm_ub_c{2} = RATE_C;
            vmap.cm_ub_a{2} = RATE_A;
            vmap.cm_ub_b{2} = RATE_B;
            vmap.cm_ub_c{2} = RATE_C;
            vmap.vad_lb{2}  = ANGMIN;
            vmap.vad_ub{2}  = ANGMAX;
            vmap.tm{2}      = TAP;
            vmap.ta{2}      = SHIFT;
            vmap.pl_fr{2}   = PF;
            vmap.ql_fr{2}   = QF;
            vmap.pl_to{2}   = PT;
            vmap.ql_to{2}   = QT;
            vmap.mu_flow_fr_ub{2} = MU_SF;
            vmap.mu_flow_to_ub{2} = MU_ST;
            vmap.mu_vad_lb{2}     = MU_ANGMIN;
            vmap.mu_vad_ub{2}     = MU_ANGMAX;
        end

        function dt = default_export_data_table(obj, spec)
            %

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, MU_ANGMAX);
        end
    end     %% methods
end         %% classdef
