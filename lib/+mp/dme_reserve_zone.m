classdef dme_reserve_zone < mp.dm_element & mp.dme_shared_opf
% mp.dme_reserve_zone - Data model element for reserve zone.
%
% Implements the data element model for reserve zone elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   =========  =========  ==================================================
%   Name       Type       Description
%   =========  =========  ==================================================
%   ``req``    *double*   zonal reserve requirement *(MW)*
%   ``zones``  *integer*  matrix defining generators included in the zone
%   ``prc``    *double*   zonal reserve price *(u/MW)* [#]_
%   =========  =========  ==================================================
%
% .. [#] Here *u* denotes the units of the objective function, e.g. USD.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        zones   % zone map for online zones / gens
        req     % reserve requirement in p.u. for each active zone
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'reserve_zone';
        end

        function label = label(obj)
            %
            label = 'Reserve Zone';
        end

        function label = labels(obj)
            %
            label = 'Reserve Zones';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'req', 'zones', 'prc'} );
        end

        function vars = export_vars(obj)
            %
            vars = {'prc'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.prc = 0;
        end

        function obj = update_status(obj, dm)
            %

            %% get gen status info
            gen_dme = dm.elements.gen;
            ng = gen_dme.nr;
            gs = gen_dme.tab.status;    %% gen status

            %% zero out columns of zones corresponding to offline gens
            z = obj.tab.zones * spdiags(gs, 0, ng, ng);

            %% update status of reserve gens
            obj.tab.status = obj.tab.status .* any(z, 2);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %
            gen_dme = dm.elements.gen;

            %% zone map and requirements for online zones / gens
            obj.zones = obj.tab.zones(obj.on, gen_dme.on);
            obj.req = obj.tab.req(obj.on) / dm.base_mva;
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                   No. of  Requirement   Price', ...
                    ' Zone ID   Status   Gens       (MW)     ($/MW)', ...
                    '---------  ------  ------  -----------  --------' } ];
            %%       12345678 -----1 12345678 1234567890.2 123456.89
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            zones = dm.elements.reserve_zone.tab.zones;
            ngr = sum(zones(k, :));
            str = sprintf('%8d %6d %8d %12.1f %9.2f', ...
                obj.tab.uid(k), obj.tab.status(k), ngr, ...
                obj.tab.req(k), obj.tab.prc(k));
        end
    end     %% methods
end         %% classdef
