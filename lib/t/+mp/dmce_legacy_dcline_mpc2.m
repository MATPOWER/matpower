classdef dmce_legacy_dcline_mpc2 < mp.dmc_element % & mp.dmce_legacy_dcline
% mp.dmce_legacy_dcline_mpc2 - Data model converter element for legacy DC line for |MATPOWER| case v2.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function name = name(obj)
            %
            name = 'legacy_dcline';
        end

        function df = data_field(obj)
            %
            df = 'dcline';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            c = idx_dcline;

            cip_fcn = @(ob, mpc, spec, vn)dcline_cost_import(ob, mpc, spec, vn);
            cep_fcn = @(ob, dme, mpc, spec, vn, ridx)dcline_cost_export(ob, dme, mpc, spec, vn, ridx);

            %% mapping for each name, default is {'col', []}
            vmap.uid                = {'IDs'};      %% consecutive IDs, starting at 1
            vmap.name               = {'cell', ''};     %% empty char
            vmap.status{2}          = c.BR_STATUS;
            vmap.source_uid         = {'cell', ''};     %% empty char
            vmap.bus_fr{2}          = c.F_BUS;
            vmap.bus_to{2}          = c.T_BUS;
            vmap.loss0{2}           = c.LOSS0;
            vmap.loss1{2}           = c.LOSS1;
            vmap.vm_setpoint_fr{2}  = c.VF;
            vmap.vm_setpoint_to{2}  = c.VT;
            vmap.p_fr_lb{2}         = c.PMIN;
            vmap.p_fr_ub{2}         = c.PMAX;
            vmap.q_fr_lb{2}         = c.QMINF;
            vmap.q_fr_ub{2}         = c.QMAXF;
            vmap.q_to_lb{2}         = c.QMINT;
            vmap.q_to_ub{2}         = c.QMAXT;
            vmap.p_fr{2}            = c.PF;
            vmap.q_fr{2}            = c.QF;
            vmap.p_to{2}            = c.PT;
            vmap.q_to{2}            = c.QT;
            if isfield(vmap, 'cost')
                vmap.cost        = {'fcn', cip_fcn, cep_fcn};
                vmap.mu_p_fr_lb{2}  = c.MU_PMIN;
                vmap.mu_p_fr_ub{2}  = c.MU_PMAX;
                vmap.mu_q_fr_lb{2}  = c.MU_QMINF;
                vmap.mu_q_fr_ub{2}  = c.MU_QMAXF;
                vmap.mu_q_to_lb{2}  = c.MU_QMINT;
                vmap.mu_q_to_ub{2}  = c.MU_QMAXT;
            end
        end

        function dt = default_export_data_table(obj, spec)
            %

            %% define named indices into data matrices
            c = idx_dcline;

            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, c.QMAXT);
        end

        function val = dcline_cost_import(obj, mpc, spec, vn)
            %
            if isfield(mpc, 'dclinecost') && spec.nr
                val = mp.dmce_gen_mpc2.gencost2cost_table(mpc.dclinecost);
            else
                val = [];
            end
        end

        function mpc = dcline_cost_export(obj, dme, mpc, spec, vn, ridx)
            %

            if dme.have_cost()
                cost = mp.dmce_gen_mpc2.cost_table2gencost( ...
                            [], dme.tab.cost, ridx);
                mpc.dclinecost(1:spec.nr, 1:size(cost, 2)) = cost;
            end
        end
    end     %% methods
end         %% classdef
