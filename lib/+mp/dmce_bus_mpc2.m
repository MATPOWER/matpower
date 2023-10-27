classdef dmce_bus_mpc2 < mp.dmc_element % & mp.dmce_bus
% mp.dmce_bus_mpc2 - Data model converter element for bus for |MATPOWER| case v2.

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
            name = 'bus';
        end

        function df = data_field(obj)
            %
            df = 'bus';
        end

        function vmap = table_var_map(obj, dme, mpc)
            %
            vmap = table_var_map@mp.dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            bni_fcn = @(ob, mpc, spec, vn)bus_name_import(ob, mpc, spec, vn, 1);
            bne_fcn = @(ob, dme, mpc, spec, vn, ridx)bus_name_export(ob, dme, mpc, spec, vn, ridx, 1);
            bsi_fcn = @(ob, mpc, spec, vn)bus_status_import(ob, mpc, spec, vn, BUS_TYPE);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}         = BUS_I;
            vmap.name           = {'fcn', bni_fcn, bne_fcn};    %% fcns to import/export from/to mpc.bus_name
            vmap.status         = {'fcn', bsi_fcn}; %% fcn w/logic for mpc.bus_types
            vmap.source_uid     = {'cell', ''};     %% empty char
            vmap.base_kv{2}     = BASE_KV;
            vmap.type{2}        = BUS_TYPE;
            vmap.area{2}        = BUS_AREA;
            vmap.zone{2}        = ZONE;
            vmap.vm_lb{2}       = VMIN;
            vmap.vm_ub{2}       = VMAX;
            vmap.va{2}          = VA;
            vmap.vm{2}          = VM;
            if isfield(vmap, 'lam_p')
                vmap.lam_p{2}       = LAM_P;
                vmap.lam_q{2}       = LAM_Q;
                vmap.mu_vm_lb{2}    = MU_VMIN;
                vmap.mu_vm_ub{2}    = MU_VMAX;
            end
        end

        function d = init_export_data(obj, dme, d, spec)
            %
            d = init_export_data@mp.dmc_element(obj, dme, d, spec); %% call parent
            if ~all(cellfun(@isempty, dme.tab.name))
                d.bus_name = cell(spec.nr, 1);
                [d.bus_name{:}] = deal('');
            end
        end

        function dt = default_export_data_table(obj, spec)
            %

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            nr = obj.default_export_data_nrows(spec);
            dt = zeros(nr, MU_VMIN);
        end

        function vals = bus_name_import(obj, mpc, spec, vn, c)
            %
            if isfield(mpc, 'bus_name')
                if spec.nr && isempty(spec.r)
                    vals = mpc.bus_name(:, c);
                else
                    vals = mpc.bus_name(spec.r, c);
                end
            else
                vals = cell(spec.nr, 1);
                [vals{:}] = deal('');
            end
        end

        function mpc = bus_name_export(obj, dme, mpc, spec, vn, ridx, c)
            %
            if isempty(ridx)
                bus_names = dme.tab.name;
            else
                bus_names = dme.tab.name(ridx, :);
            end
            if ~all(cellfun(@isempty, bus_names))
                if spec.nr && isempty(spec.r)
                    if isempty(ridx)
                        mpc.bus_name(:, c) = bus_names;
                    else
                        mpc.bus_name(ridx, c) = bus_names;
                    end
                else
                    if isempty(ridx)
                        mpc.bus_name(spec.r, c) = bus_names;
                    else
                        mpc.bus_name(spec.r(ridx), c) = bus_names;
                    end
                end
            end
        end

        function vals = bus_status_import(obj, mpc, spec, vn, c)
            %

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

            if spec.nr && isempty(spec.r)
                vals = mpc.bus(:, BUS_TYPE) ~= NONE;
            else
                vals = mpc.bus(spec.r, BUS_TYPE) ~= NONE;
            end
        end
    end     %% methods
end         %% classdef
