classdef dme_load3p < mp.dm_element
% mp.dme_load3p - Data model element for 3-phase load.
%
% Implements the data element model for 3-phase load elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ================  =========  =====================================
%   Name              Type       Description
%   ================  =========  =====================================
%   ``bus``           *integer*  bus ID (``uid``) of 3-phase bus
%   ``pd1``           *double*   phase 1 active power demand *(kW)*
%   ``pd2``           *double*   phase 2 active power demand *(kW)*
%   ``pd3``           *double*   phase 3 active power demand *(kW)*
%   ``pf1``           *double*   phase 1 power factor
%   ``pf2``           *double*   phase 2 power factor
%   ``pf3``           *double*   phase 3 power factor
%   ================  =========  =====================================

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     % bus index vector (all loads)
        pd1     % phase 1 active power demand (p.u.) for loads that are on
        pd2     % phase 2 active power demand (p.u.) for loads that are on
        pd3     % phase 3 active power demand (p.u.) for loads that are on
        pf1     % phase 1 power factor for loads that are on
        pf2     % phase 2 power factor for loads that are on
        pf3     % phase 3 power factor for loads that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'load3p';
        end

        function label = label(obj)
            %
            label = '3-ph Load';
        end

        function label = labels(obj)
            %
            label = '3-ph Loads';
        end

        function name = cxn_type(obj)
            %
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            %
            name = 'bus';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'pd1', 'pd2', 'pd3', 'pf1', 'pf2', 'pf3'});
        end

%         function vars = export_vars(obj)
%             vars = {};
%         end

%         function s = export_vars_offline_val(obj)
%             s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
%         end

        function obj = initialize(obj, dm)
            %
            initialize@mp.dm_element(obj, dm);  %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus3p.ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bs = dm.elements.bus3p.tab.status;  %% bus status

            %% update status of loads at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %
            obj.pd1 = obj.tab.pd1(obj.on) / dm.base_kva;
            obj.pd2 = obj.tab.pd2(obj.on) / dm.base_kva;
            obj.pd3 = obj.tab.pd3(obj.on) / dm.base_kva;
            obj.pf1 = obj.tab.pf1(obj.on);
            obj.pf2 = obj.tab.pf2(obj.on);
            obj.pf3 = obj.tab.pf3(obj.on);
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print generation summary
            pd = [ obj.tab.pd1 obj.tab.pd2 obj.tab.pd3 ];
            pf = [ obj.tab.pf1 obj.tab.pf2 obj.tab.pf3 ];
            qd = pd .* tan(acos(pf));
            fprintf(fd, '  %-29s %12.1f kW %12.1f kVAr\n', 'Total 3-ph load', ...
                sum(sum(pd(obj.on, :))), sum(sum(qd(obj.on, :))) );
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '  3-ph      3-ph             Phase A Power     Phase B Power     Phase C Power', ...
                    'Load ID    Bus ID   Status   (kW)     (PF)     (kW)     (PF)     (kW)     (PF)', ...
                    '--------  --------  ------  -------  ------   -------  ------   -------  ------' } ];
            %%       1234567 123456789 -----1 1234567.90 12.4567 123456.89 12.4567 123456.89 12.4567
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            str = sprintf('%7d %9d %6d %10.2f %7.4f %9.2f %7.4f %9.2f %7.4f', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), ...
                obj.tab.pd1(k), obj.tab.pf1(k), ...
                obj.tab.pd2(k), obj.tab.pf2(k), ...
                obj.tab.pd3(k), obj.tab.pf3(k) );
        end
    end     %% methods
end         %% classdef
