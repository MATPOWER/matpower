classdef dme_xfmr3p < mp.dm_element
% mp.dme_xfmr3p - Data model element for 3-phase transformer.
%
% Implements the data element model for 3-phase transformer elements.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ============  =========  =============================================
%   Name          Type       Description
%   ============  =========  =============================================
%   ``bus_fr``    *integer*  bus ID (``uid``) of "from" 3-phase bus
%   ``bus_to``    *integer*  bus ID (``uid``) of "to" 3-phase bus
%   ``r``         *double*   series resistance *(p.u.)*
%   ``x``         *double*   series reactance *(p.u.)*
%   ``base_kva``  *double*   transformer kVA base *(kVA)*
%   ``base_kv``   *double*   transformer kV base *(kV)*
%   ``pl1_fr``    *double*   phase 1 active power injection at "from" end *(kW)*
%   ``ql1_fr``    *double*   phase 1 reactive power injection at "from" end *(kVAr)*
%   ``pl2_fr``    *double*   phase 2 active power injection at "from" end *(kW)*
%   ``ql2_fr``    *double*   phase 2 reactive power injection at "from" end *(kVAr)*
%   ``pl3_fr``    *double*   phase 3 active power injection at "from" end *(kW)*
%   ``ql3_fr``    *double*   phase 3 reactive power injection at "from" end *(kVAr)*
%   ``pl1_to``    *double*   phase 1 active power injection at "to" end *(kW)*
%   ``ql1_to``    *double*   phase 1 reactive power injection at "to" end *(kVAr)*
%   ``pl2_to``    *double*   phase 2 active power injection at "to" end *(kW)*
%   ``ql2_to``    *double*   phase 2 reactive power injection at "to" end *(kVAr)*
%   ``pl3_to``    *double*   phase 3 active power injection at "to" end *(kW)*
%   ``ql3_to``    *double*   phase 3 reactive power injection at "to" end *(kVAr)*
%   ============  =========  =============================================

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    % bus index vector for "from" bus (all transformers)
        tbus    % bus index vector for "to" bus (all transformers)
        r       % series resistance (p.u.) for transformers that are on
        x       % series reactance (p.u.) for transformers that are on
        base_kva%
        base_kv %
%         g_to    %% shunt conductance (p.u.) at "to" end for transformers that are on
%         b_fr    %% shunt susceptance (p.u.) at "from" end for transformers that are on
%         b_to    %% shunt susceptance (p.u.) at "to" end for transformers that are on
%         tm      %% transformer off-nominal turns ratio for transformers that are on
%         ta      %% xformer phase-shift angle (radians) for transformers that are on
%         rate_a  %% long term flow limit (p.u.) for transformers that are on
    end     %% properties

    methods
        function name = name(obj)
            %
            name = 'xfmr3p';
        end

        function label = label(obj)
            %
            label = '3-ph Transformer';
        end

        function label = labels(obj)
            %
            label = '3-ph Transformers';
        end

        function name = cxn_type(obj)
            %
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            %
            name = {'fbus', 'tbus'};
        end

        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'base_kva', 'base_kv', ...
                 'pl1_fr', 'ql1_fr', 'pl2_fr', 'ql2_fr', 'pl3_fr', 'ql3_fr', ...
                 'pl1_to', 'ql1_to', 'pl2_to', 'ql2_to', 'pl3_to', 'ql3_to' ...
                 });
        end

%         function vars = export_vars(obj)
%             vars = {'pl1_fr', 'ql1_fr', ...
%                     'pl2_fr', 'ql2_fr', ...
%                     'pl3_fr', 'ql3_fr', ...
%                     'pl1_to', 'ql1_to', ...
%                     'pl2_to', 'ql2_to', ...
%                     'pl3_to', 'ql3_to' };
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
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %

            %% get bus status info
            bs = dm.elements.bus3p.tab.status;  %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.fbus) & ...
                                              bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %

            obj.r = obj.tab.r(obj.on);
            obj.x = obj.tab.x(obj.on);
            obj.base_kva = obj.tab.base_kva(obj.on);
            obj.base_kv  = obj.tab.base_kv(obj.on);
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
            %
            switch section
                case 'det'
                    %% compute currents/powers for pp_args
                    s_fr = [  obj.tab.pl1_fr + 1j * obj.tab.ql1_fr ...
                            obj.tab.pl2_fr + 1j * obj.tab.ql2_fr ...
                            obj.tab.pl3_fr + 1j * obj.tab.ql3_fr    ];
                    s_to = [  obj.tab.pl1_to + 1j * obj.tab.ql1_to ...
                            obj.tab.pl2_to + 1j * obj.tab.ql2_to ...
                            obj.tab.pl3_to + 1j * obj.tab.ql3_to    ];
                    t = dm.elements.bus3p.tab;    %% bus3p table
                    vm = [ t.vm1 t.vm2 t.vm3 ] .* (t.base_kv/sqrt(3) * [1 1 1]);
                    va = [ t.va1 t.va2 t.va3 ];
                    v_ = vm .* exp(1j * va * pi/180);
                    i_fr = conj( s_fr ./ v_(obj.fbus, :));
                    i_to = conj( s_to ./ v_(obj.tbus, :));
                    c_fr = struct('cm', abs(i_fr), 'ca', angle(i_fr) * 180/pi);
                    c_to = struct('cm', abs(i_to), 'ca', angle(i_to) * 180/pi);

                    pp_args.xfmr3p = {'c', 'f', abs(i_fr), angle(i_fr) * 180/pi};
                    pretty_print@mp.dm_element(obj, dm, section, out_e, mpopt, fd, pp_args);

                    pp_args.xfmr3p = {'c', 't', abs(i_to), angle(i_to) * 180/pi};
                    pretty_print@mp.dm_element(obj, dm, section, out_e, mpopt, fd, pp_args);

                    pp_args.xfmr3p = {'s', 'f', real(s_fr), imag(s_fr)};
                    pretty_print@mp.dm_element(obj, dm, section, out_e, mpopt, fd, pp_args);

                    pp_args.xfmr3p = {'s', 't', real(s_to), imag(s_to)};
                    pretty_print@mp.dm_element(obj, dm, section, out_e, mpopt, fd, pp_args);
                otherwise
                    pretty_print@mp.dm_element(obj, dm, section, out_e, mpopt, fd, pp_args);
            end
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
            t = obj.tab;
            ploss = [ t.pl1_fr t.pl2_fr t.pl3_fr ] + ...
                    [ t.pl1_to t.pl2_to t.pl3_to ];
            qloss = [ t.ql1_fr t.ql2_fr t.ql3_fr ] + ...
                    [ t.ql1_to t.ql2_to t.ql3_to ];
            fprintf(fd, '  %-29s %12.1f kW %12.1f kVAr\n', 'Total 3-ph transformer loss', ...
                sum(sum(ploss(obj.on, :))), sum(sum(qloss(obj.on, :))) );
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            cs = pp_args.xfmr3p{1};
            ft = pp_args.xfmr3p{2};
            if cs == 'c' && ft == 'f'
                h1 = pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args);
            else
                h1 = {};
            end
            if cs == 'c'
                if ft == 'f'
                    h2 = {'-->  Current Injections at "From" Bus' };
                else    %% ft == 't'
                    h2 = {'', ...
                        '<--  Current Injections at "To" Bus' };
                end
            else    %% cs == 's'
                if ft == 'f'
                    h2 = {'', ...
                        '-->  Power Injections at "From" Bus' };
                else    %% ft == 't'
                    h2 = {'', ...
                        '<--  Power Injections at "To" Bus' };
                end
            end
            switch cs
                case 'c'
                    h = [ h1 h2 ...
                        {   '  3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current', ...
                            'Xfrm ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' } ];
                    %%       1234567 123456789 123456789 -----1 123456.89 12345.7 12345.78 12345.7 12345.78 12345.7
                case 's'
                    h = [ h1 h2 ...
                        {   '  3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power', ...
                            'Xfmr ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' } ];
                    %%       1234567 123456789 123456789 -----1 1234567.9 12345.7 123456.8 12345.7 123456.8 12345.7
            end
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            switch pp_args.xfmr3p{1}    %% cs
                case 'c'
                    cm = pp_args.xfmr3p{3}(k, :);
                    ca = pp_args.xfmr3p{4}(k, :);
                    str = sprintf('%7d %9d %9d %6d %9.2f %7.1f %8.2f %7.1f %8.2f %7.1f', ...
                        obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                        obj.tab.status(k), ...
                        cm(:, 1), ca(:, 1), ...
                        cm(:, 2), ca(:, 2), ...
                        cm(:, 3), ca(:, 3) );
                case 's'
                    p = pp_args.xfmr3p{3}(k, :);
                    q = pp_args.xfmr3p{4}(k, :);
                    str = sprintf('%7d %9d %9d %6d %9.1f %7.1f %8.1f %7.1f %8.1f %7.1f', ...
                        obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                        obj.tab.status(k), ...
                        p(:, 1), q(:, 1), ...
                        p(:, 2), q(:, 2), ...
                        p(:, 3), q(:, 3) );
            end
        end
    end     %% methods
end         %% classdef
