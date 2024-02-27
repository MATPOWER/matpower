classdef dme_bus_opf < mp.dme_bus & mp.dme_shared_opf
% mp.dme_bus_opf - Data model element for bus for OPF.
%
% To parent class :class:`mp.dme_bus`, adds shadow prices on power balance
% and voltage magnitude limits, and pretty-printing for **lim** sections.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   ============  ========  ==============================================
%   Name          Type      Description
%   ============  ========  ==============================================
%   ``lam_p``     *double*  active power nodal price, i.e. shadow price on
%                           active power balance constraint *(u/MW)* [#]_
%   ``lam_q``     *double*  reactive power nodal price, i.e. shadow price on
%                           reactive power balance constraint *(u/MVAr)* [1]_
%   ``mu_vm_lb``  *double*  shadow price on voltage magnitude lower
%                           bound *(u/p.u.)* [1]_
%   ``mu_vm_ub``  *double*  shadow price on voltage magnitude upper
%                           bound *(u/p.u.)* [1]_
%   ============  ========  ==============================================
%
% .. [#] Here *u* denotes the units of the objective function, e.g. USD.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function names = main_table_var_names(obj)
            %
            names = horzcat( main_table_var_names@mp.dme_bus(obj), ...
                {'lam_p', 'lam_q', 'mu_vm_lb', 'mu_vm_ub'});
        end

        function vars = export_vars(obj)
            %
            vars = horzcat( export_vars@mp.dme_bus(obj), ...
                {'vm_lb', 'vm_ub', 'lam_p', 'lam_q', ...
                    'mu_vm_lb', 'mu_vm_ub'} );
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dme_bus(obj);    %% call parent
            s.lam_p = 0;
            s.lam_q = 0;
            s.mu_vm_lb = 0;
            s.mu_vm_ub = 0;
        end

        function obj = pp_data_ext(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %

            %% call parent
            pp_data_ext@mp.dme_bus(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print bus extremes
            [min_lam_p, min_lam_p_i] = min(obj.tab.lam_p);
            [max_lam_p, max_lam_p_i] = max(obj.tab.lam_p);

            fprintf(fd, '  %-29s %15s @ %11s %16s @ %s\n', ...
                'Lambda P (LMP active power)', ...
                sprintf('%8.2f $/MWh', min_lam_p), ...
                sprintf('bus %-7d', obj.tab.uid(min_lam_p_i)), ...
                sprintf('%8.2f $/MWh', max_lam_p), ...
                sprintf('bus %d', obj.tab.uid(max_lam_p_i)) );

            if mpopt.model(1) ~= 'D'    %% AC model
                [min_lam_q, min_lam_q_i] = min(obj.tab.lam_q);
                [max_lam_q, max_lam_q_i] = max(obj.tab.lam_q);
                fprintf(fd, '  %-29s %15s @ %11s %16s @ %s\n', ...
                    'Lambda Q (LMP reactive power)', ...
                    sprintf('%7.2f $/MVArh', min_lam_q), ...
                    sprintf('bus %-7d', obj.tab.uid(min_lam_q_i)), ...
                    sprintf('%7.2f $/MVArh', max_lam_q), ...
                    sprintf('bus %d', obj.tab.uid(max_lam_q_i)) );
            end
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            %
            if mpopt.model(1) ~= 'D'    %% AC model
                h = pp_get_headers_det@mp.dme_bus(obj, dm, out_e, mpopt, pp_args);
                h{end-2} = [ h{end-2} '            Lambda (LMP)'];
                h{end-1} = [ h{end-1} '  P($/MWh)  Q($/MVAr-hr)'];
                h{end}   = [ h{end}   '  --------  ------------'];
            else                        %% DC model
                h = pp_get_headers_det@mp.dme_bus(obj, dm, out_e, mpopt, pp_args);
                h{end-2} = [ h{end-2} '       Lambda (LMP)'];
                h{end-1} = [ h{end-1} '  P($/MWh)'];
                h{end}   = [ h{end}   '  --------'];
            end
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if mpopt.model(1) ~= 'D'    %% AC model
                str = [ ...
                    pp_data_row_det@mp.dme_bus(obj, dm, k, out_e, mpopt, fd, pp_args) ...
                    sprintf(' %9.3f %13.3f', ...
                        obj.tab.lam_p(k), obj.tab.lam_q(k))];
            else                        %% DC model
                str = [ ...
                    pp_data_row_det@mp.dme_bus(obj, dm, k, out_e, mpopt, fd, pp_args) ...
                    sprintf(' %9.3f', obj.tab.lam_p(k))];
            end
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            rows = find( obj.tab.status & ( ...
                        obj.tab.vm < obj.tab.vm_lb + obj.ctol | ...
                        obj.tab.vm > obj.tab.vm_ub - obj.ctol | ...
                        obj.tab.mu_vm_lb > obj.ptol | ...
                        obj.tab.mu_vm_ub > obj.ptol ));
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            h = [ pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args) ...
                {   '                     Voltage Magnitude Limits', ...
                    ' Bus ID     mu LB      LB       vm       UB       mu UB', ...
                    '--------  ---------  -------  -------  -------   --------' } ];
            %%       1234567  12345.789    1.345    1.345    1.345  12345.789
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            if obj.tab.vm(k) < obj.tab.vm_lb(k) + obj.ctol || ...
                    obj.tab.mu_vm_lb(k) > obj.ptol
                mu_lb = sprintf('%10.3f', obj.tab.mu_vm_lb(k));
            else
                mu_lb = '      -   ';
            end
            if obj.tab.vm(k) > obj.tab.vm_ub(k) - obj.ctol || ...
                    obj.tab.mu_vm_ub(k) > obj.ptol
                mu_ub = sprintf('%10.3f', obj.tab.mu_vm_ub(k));
            else
                mu_ub = '      -   ';
            end

            str = sprintf('%7d %10s %8.3f %8.3f %8.3f %10s', ...
                obj.tab.uid(k), mu_lb, obj.tab.vm_lb(k), ...
                obj.tab.vm(k), obj.tab.vm_ub(k), mu_ub);
        end
    end     %% methods
end         %% classdef
