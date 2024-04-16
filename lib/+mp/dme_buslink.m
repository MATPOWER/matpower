classdef dme_buslink < mp.dm_element
% mp.dme_buslink - Data model element for 1-to-3-phase buslink.
%
% Implements the data element model for 1-to-3-phase buslink elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ===========  =========  ========================================
%   Name         Type       Description
%   ===========  =========  ========================================
%   ``bus``      *integer*  bus ID (``uid``) of single phase bus
%   ``bus3p``    *integer*  bus ID (``uid``) of 3-phase bus
%   ===========  =========  ========================================

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         % bus index vector (all buslinks)
        bus3p       % bus3p index vector (all buslinks)
        pg1_start   % initial phase 1 active power (p.u.) for buslinks that are on
        pg2_start   % initial phase 2 active power (p.u.) for buslinks that are on
        pg3_start   % initial phase 3 active power (p.u.) for buslinks that are on
        qg1_start   % initial phase 1 reactive power (p.u.) for buslinks that are on
        qg2_start   % initial phase 2 reactive power (p.u.) for buslinks that are on
        qg3_start   % initial phase 3 reactive power (p.u.) for buslinks that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'buslink';
        end

        function label = label(obj)
            %
            label = 'Bus Link';
        end

        function label = labels(obj)
            %
            label = 'Bus Links';
        end

        function name = cxn_type(obj)
            %
            name = {'bus', 'bus3p'};
        end

        function name = cxn_idx_prop(obj)
            %
            name = {'bus', 'bus3p'};
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'bus3p'});
        end

%         function vars = export_vars(obj)
%             vars = {};
%         end

%         function s = export_vars_offline_val(obj)
%             s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
%         end

        function obj = initialize(obj, dm)
            %
            initialize@mp.dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i  = dm.elements.bus.ID2i;    %% bus num to idx mapping
            b2i3 = dm.elements.bus3p.ID2i;  %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus   = b2i( obj.tab.bus);
            obj.bus3p = b2i3(obj.tab.bus3p);
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bs  = dm.elements.bus.tab.status;   %% bus status
            bs3 = dm.elements.bus3p.tab.status; %% bus3p status

            %% update status of buslinks at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus) & bs3(obj.bus3p);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %

            %% check for matching base_kv
            base_kv1 = dm.elements.bus.tab.base_kv(obj.bus);
            base_kv3 = dm.elements.bus3p.tab.base_kv(obj.bus3p);
            if any(base_kv1 ~= base_kv3)
                error('mp.dme_buslink.build_params: buslink objects can only link buses with identical base_kv');
            end

            %% check for at least one PQ bus per buslink
            type1 = dm.elements.bus.tab.type(obj.bus);
            type3 = dm.elements.bus3p.tab.type(obj.bus3p);
            if any(type1 ~= mp.NODE_TYPE.PQ & type3 ~= mp.NODE_TYPE.PQ)
                error('mp.dme_buslink.build_params: at lease one of the buses linked by a buslink must be of type PQ');
            end

            %% check for balanced voltage magnitudes for REF and PV buses
            %% on 3-phase side
            refpv = (type3 == mp.NODE_TYPE.REF | type3 == mp.NODE_TYPE.PV);
            vm1_ref_pv = dm.elements.bus3p.tab.vm1(obj.bus3p(refpv));
            vm2_ref_pv = dm.elements.bus3p.tab.vm2(obj.bus3p(refpv));
            vm3_ref_pv = dm.elements.bus3p.tab.vm3(obj.bus3p(refpv));
            if any(vm1_ref_pv ~= vm2_ref_pv | vm2_ref_pv ~= vm3_ref_pv)
                error('mp.dme_buslink.build_params: buslink objects can only link to REF or PV buses with balanced voltage magnitudes');
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                      3-ph', ...
                    'Link ID    Bus ID    Bus ID   Status', ...
                    '--------  --------  --------  ------' } ];
            %%       1234567 123456789 123456789 -----1
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            str = sprintf('%7d %9d %9d %6d', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.bus3p(k), obj.tab.status(k) );
        end
    end     %% methods
end         %% classdef
