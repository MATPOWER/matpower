classdef dm_converter_mpc2 < mp.dm_converter
% mp.dm_converter_mpc2 - |MATPOWER| **data model converter** for MATPOWER case v2.
%
% This class implements importing/exporting of data models for version 2
% of the classic |MATPOWER| case format. That is, the *data source* **d**
% for this class is expected to be a |MATPOWER| case struct.
%
% mp.dm_converter_mpc2 Methods:
%   * dm_converter_mpc2 - constructor
%   * format_tag - return char array identifier for data source/format (``'mpc2'``)
%   * import - import data from a |MATPOWER| case struct into a data model
%   * export - export data from a data model to a |MATPOWER| case struct
%   * save - save |MATPOWER| case struct to a file
%
% See also mp.dm_converter.

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
        function obj = dm_converter_mpc2()
            % Specify the element classes for handling |MATPOWER| case format.

            %% call parent constructor
            obj@mp.dm_converter();
            obj.element_classes = ...
                { @mp.dmce_bus_mpc2, @mp.dmce_gen_mpc2, @mp.dmce_load_mpc2, ...
                    @mp.dmce_branch_mpc2, @mp.dmce_shunt_mpc2 };
        end

        function tag = format_tag(obj)
            % Return identifier tag ``'mpc2'`` for version 2 |MATPOWER| case format.

            tag = 'mpc2';
        end

        function dm = import(obj, dm, d)
            % Import data from a version 2 |MATPOWER| case struct into a data model.

            if ~isstruct(d)
                d = loadcase(d);
            end
            if isfield(d, 'baseMVA');
                dm.base_mva = d.baseMVA;
            elseif isfield(d, 'bus') && ~isempty(d.bus)
                error('mp.dm_converter_mpc2.import: ''baseMVA'' must be defined for a case with a ''bus'' field.');
            end
            if isfield(d, 'basekVA');
                dm.base_kva = d.basekVA;
            elseif isfield(d, 'bus3p') && ~isempty(d.bus3p)
                error('mp.dm_converter_mpc2.import: ''basekVA'' must be defined for a case with a ''bus3p'' field.');
            end
            dm.source = d;
            dm = import@mp.dm_converter(obj, dm, d);
        end

        function d = init_export(obj, dm)
            % Initialize a |MATPOWER| case struct for export.

            d = struct('version', '2', 'baseMVA', dm.base_mva);
        end

        function fname_out = save(obj, fname, d)
            % Save a |MATPOWER| case struct to a file.

            fname_out = savecase(fname, d);
        end
    end     %% methods
end         %% classdef
