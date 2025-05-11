classdef dme_shunt3p < mp.dm_element
% mp.dme_shunt3p - Data model element for 3-phase shunt.
%
% Implements the data element model for 3-phase shunt elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ========  =========  =====================================
%   Name      Type       Description
%   ========  =========  =====================================
%   ``bus``   *integer*  bus ID (``uid``)
%   ``gs1``   *double*   phase 1 shunt conductance, specified as
%                        nominal [#]_ active power demand *(MW)*
%   ``gs2``   *double*   phase 2 shunt conductance, specified as
%                        nominal [#]_ active power demand *(MW)*
%   ``gs3``   *double*   phase 3 shunt conductance, specified as
%                        nominal [#]_ active power demand *(MW)*
%   ``bs1``   *double*   phase 1 shunt susceptance, specified as
%                        nominal [1]_ reactive power injection *(MVAr)*
%   ``bs2``   *double*   phase 1 shunt susceptance, specified as
%                        nominal [1]_ reactive power injection *(MVAr)*
%   ``bs3``   *double*   phase 3 shunt susceptance, specified as
%                        nominal [1]_ reactive power injection *(MVAr)*
%   ========  =========  =====================================
%
% .. [#] *Nominal* means for a voltage of 1 p.u.

%   MATPOWER
%   Copyright (c) 2021-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     % bus index vector (all loads)
        gs1     % phase 1 shunt conductance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        gs2     % phase 2 shunt conductance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        gs3     % phase 3 shunt conductance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        bs1     % phase 1 shunt susceptance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        bs2     % phase 2 shunt susceptance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
        bs3     % phase 3 shunt susceptance (p.u. active power demanded at
                % V = 1.0 p.u.) for shunts that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'shunt3p';
        end

        function label = label(obj)
            %
            label = '3-ph Shunt';
        end

        function label = labels(obj)
            %
            label = '3-ph Shunts';
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
                {'bus', 'gs1', 'gs2', 'gs3', 'bs1', 'bs2', 'bs3', ...
                 'p1', 'p2', 'p3', 'q1', 'q2', 'q3'});
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
            obj.gs1 = obj.tab.gs1(obj.on) / dm.base_kva;
            obj.gs2 = obj.tab.gs2(obj.on) / dm.base_kva;
            obj.gs3 = obj.tab.gs3(obj.on) / dm.base_kva;
            obj.bs1 = obj.tab.bs1(obj.on) / dm.base_kva;
            obj.bs2 = obj.tab.bs2(obj.on) / dm.base_kva;
            obj.bs3 = obj.tab.bs3(obj.on) / dm.base_kva;
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
            p = [ obj.tab.p1 obj.tab.p2 obj.tab.p3 ];
            q = [ obj.tab.q1 obj.tab.q2 obj.tab.q3 ];

            fprintf(fd, '  %-29s %12.1f kW', 'Total shunt', sum(sum(p)));
            if mpopt.model(1) ~= 'D'    %% AC model
                fprintf(fd, ' %12.1f kVAr', sum(sum(q)));
            end
            fprintf(fd, '\n');
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f kW', '  online', ...
                                                sum(sum(p(obj.on, :))));
                if mpopt.model(1) ~= 'D'    %% AC model
                    fprintf(fd, ' %12.1f kVAr', sum(sum(q(obj.on, :))));
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
                {   '                                                                Power consumption', ...
                    '   3-ph      3-ph             Phase A Power      Phase B Power      Phase C Power', ...
                    'Shunt ID    Bus ID   Status   (kW)    (kVAr)     (kW)    (kVAr)     (kW)     (kVAr)', ...
                    '---------  --------  ------  -------  -------   -------  -------   -------  -------' } ];
            %%       12345678 123456789 -----1 123456.89 123456.89 123456.89 123456.89 123456.89 123456.89
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if obj.tab.status(k) && abs(obj.tab.p1(k)) > 1e-5
                p1 = sprintf('%9.2f', obj.tab.p1(k));
            else
                p1 = '      -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.q1(k)) > 1e-5
                q1 = sprintf('%9.2f', obj.tab.q1(k));
            else
                q1 = '      -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.p2(k)) > 1e-5
                p2 = sprintf('%9.2f', obj.tab.p2(k));
            else
                p2 = '      -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.q2(k)) > 1e-5
                q2 = sprintf('%9.2f', obj.tab.q2(k));
            else
                q2 = '      -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.p3(k)) > 1e-5
                p3 = sprintf('%9.2f', obj.tab.p3(k));
            else
                p3 = '      -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.q3(k)) > 1e-5
                q3 = sprintf('%9.2f', obj.tab.q3(k));
            else
                q3 = '      -  ';
            end
            str = sprintf('%8d %9d %6d %9s %9s %9s %9s %9s %9s', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), ...
                                           p1, q1, p2 ,q2, p3, q3);
        end        
    end     %% methods
end         %% classdef
