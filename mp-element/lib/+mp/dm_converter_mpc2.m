classdef dm_converter_mpc2 < mp.dm_converter
%MP.DM_CONVERTER_MPC2  MATPOWER data model converter for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = dm_converter_mpc2()
            %% call parent constructor
            obj@mp.dm_converter();
            obj.element_classes = ...
                { @mp.dmce_bus_mpc2, @mp.dmce_gen_mpc2, @mp.dmce_load_mpc2, ...
                    @mp.dmce_branch_mpc2, @mp.dmce_shunt_mpc2 };
        end

        function tag = format_tag(obj)
            tag = 'mpc2';
        end

        function dm = import(obj, dm, d)
            if ~isstruct(d)
                d = loadcase(d);
            end
            if isfield(d, 'baseMVA');
                dm.base_mva = d.baseMVA;
            elseif isfield(d, 'bus') && ~isempty(d.bus)
                error('mp.dm_converter_mpc2/import: ''baseMVA'' must be defined for a case with a ''bus'' field.');
            end
            if isfield(d, 'basekVA');
                dm.base_kva = d.basekVA;
            elseif isfield(d, 'bus3p') && ~isempty(d.bus3p)
                error('mp.dm_converter_mpc2/import: ''basekVA'' must be defined for a case with a ''bus3p'' field.');
            end
            dm.source = d;
            dm = import@mp.dm_converter(obj, dm, d);
        end

        function d = init_export(obj, dm)
            d = struct('version', '2', 'baseMVA', dm.base_mva);
        end

        function fname_out = save(obj, fname, d)
            fname_out = savecase(fname, d);
        end
    end     %% methods
end         %% classdef
