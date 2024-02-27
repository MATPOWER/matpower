classdef dme_shunt < mp.dm_element
% mp.dme_shunt - Data model element for shunt.
%
% Implements the data element model for shunt elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ========  =========  =====================================
%   Name      Type       Description
%   ========  =========  =====================================
%   ``bus``   *integer*  bus ID (``uid``)
%   ``gs``    *double*   :math:`g_s`, shunt conductance, specified as
%                        nominal [#]_ active power demand *(MW)*
%   ``bs``    *double*   :math:`b_s`, shunt susceptance, specified as
%                        nominal [1]_ reactive power injection *(MVAr)*
%   ``p``     *double*   :math:`p`, total active power absorbed *(MW)*
%   ``q``     *double*   :math:`q`, total reactive power absorbed *(MVAr)*
%   ========  =========  =====================================
%
% .. [#] *Nominal* means for a voltage of 1 p.u.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     % bus index vector (all shunts)
        gs      % shunt conductance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        bs      % shunt susceptance (p.u. reactive power injected at
                % V = 1.0 p.u.) for shunts that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'shunt';
        end

        function label = label(obj)
            %
            label = 'Fixed Shunt';
        end

        function label = labels(obj)
            %
            label = 'Fixed Shunts';
        end

        function name = cxn_type(obj)
            %
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            %
            name = 'bus';
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'gs', 'bs', 'p', 'q'});
        end

%         function vars = export_vars(obj)
%             %
%             vars = {};
%         end
% 
%         function s = export_vars_offline_val(obj)
%             %
%             s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
%         end

        function nr = count(obj, dm)
            %
            nr = count@mp.dm_element(obj, dm);
            if nr
                obj.bus = obj.tab.source_uid;
            end
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bs = dm.elements.bus.tab.status;    %% bus status

            %% update status of shunts at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %
            obj.gs = obj.tab.gs(obj.on) / dm.base_mva;
            obj.bs = obj.tab.bs(obj.on) / dm.base_mva;
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print shunt summary
            fprintf(fd, '  %-29s %12.1f MW', 'Total shunt', ...
                                            sum(obj.tab.p));
            if mpopt.model(1) ~= 'D'    %% AC model
                fprintf(fd, ' %12.1f MVAr', sum(obj.tab.q));
            end
            fprintf(fd, '\n');
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW', '  online', ...
                                                sum(obj.tab.p(obj.on)));
                if mpopt.model(1) ~= 'D'    %% AC model
                    fprintf(fd, ' %12.1f MVAr', sum(obj.tab.q(obj.on)));
                end
                fprintf(fd, '\n');
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                             Power Consumption', ...
                    'Shunt ID   Bus ID   Status   P (MW)   Q (MVAr)', ...
                    '--------  --------  ------  --------  --------' } ];
            %%       1234567 123456789 -----1 1234567.90 123456.89
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            %
            f = {'                            --------  --------',
                sprintf('%18s Total:%10.2f %9.2f', ...
                    '', sum(obj.tab.p(obj.on)), sum(obj.tab.q(obj.on)))};
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if obj.tab.status(k) && abs(obj.tab.p(k)) > 1e-5
                p = sprintf('%10.2f', obj.tab.p(k));
            else
                p = '       -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.q(k)) > 1e-5
                q = sprintf('%9.2f', obj.tab.q(k));
            else
                q = '      -  ';
            end
            str = sprintf('%7d %9d %6d %10s %9s', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), p, q);
        end
    end     %% methods
end         %% classdef
