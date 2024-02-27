classdef dme_branch_opf < mp.dme_branch & mp.dme_shared_opf
% mp.dme_branch_opf - Data model element for branch for OPF.
%
% To parent class :class:`mp.dme_branch`, adds shadow prices on flow and
% angle difference limits, and pretty-printing for **lim** sections.
%
% Adds the following columns in the main data table, found in the
% :attr:`tab` property:
%
%   =================  ========  ===================================
%   Name               Type      Description
%   =================  ========  ===================================
%   ``mu_flow_fr_ub``  *double*  shadow price on flow constraint at
%                                "from" end *(u/MVA)* [#]_
%   ``mu_flow_to_ub``  *double*  shadow price on flow constraint at
%                                "to" end *(u/MVA)* [1]_
%   ``mu_vad_lb``      *double*  shadow price on lower bound of
%                                voltage angle difference constraint
%                                *(u/degree)* [1]_
%   ``mu_vad_ub``      *double*  shadow price on upper bound of
%                                voltage angle difference constraint
%                                *(u/degree)* [1]_
%   =================  ========  ===================================
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
            names = horzcat( main_table_var_names@mp.dme_branch(obj), ...
                {'mu_flow_fr_ub', 'mu_flow_to_ub', ...
                 'mu_vad_lb', 'mu_vad_ub'} );
%                 'mu_sm_fr_ub', 'mu_sm_to_ub', ...
%                 'mu_pl_fr_ub', 'mu_pl_to_ub', ...
%                 'mu_cm_fr_ub', 'mu_cm_to_ub', ...
        end

        function vars = export_vars(obj)
            %
            vars = horzcat( export_vars@mp.dme_branch(obj), ...
                {'mu_flow_fr_ub', 'mu_flow_to_ub', ...
                 'mu_vad_lb', 'mu_vad_ub'} );
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dme_branch(obj);     %% call parent
            s.mu_flow_fr_ub = 0;
            s.mu_flow_to_ub = 0;
            s.mu_vad_lb = 0;
            s.mu_vad_ub = 0;
        end

        function obj = pretty_print(obj, dm, section, out_e, mpopt, fd, pp_args)
            %
            switch section
                case 'lim'
                    %% compute flows and limits to pass to parent
                    switch upper(mpopt.opf.flow_lim(1))
                        case {'P', '2'}     %% active power
                            fm_fr = obj.tab.pl_fr;
                            fm_to = obj.tab.pl_to;
                            fm_ub = obj.tab.sm_ub_a;
                        case 'I'            %% current
                            vm = dm.elements.bus.tab.vm;
                            va = dm.elements.bus.tab.va;
                            v_ = vm .* exp(1j * va * pi/180);
                            fm_fr = abs((obj.tab.pl_fr + 1j * obj.tab.ql_fr) ./ v_(obj.fbus));
                            fm_to = abs((obj.tab.pl_to + 1j * obj.tab.ql_to) ./ v_(obj.tbus));
                            fm_ub = obj.tab.cm_ub_a;
                        otherwise           %% apparent power
                            fm_fr = sqrt(obj.tab.pl_fr.^2 + obj.tab.ql_fr.^2);
                            fm_to = sqrt(obj.tab.pl_to.^2 + obj.tab.ql_to.^2);
                            fm_ub = obj.tab.sm_ub_a;
                    end
                    pp_args.branch.flow_magnitude = struct( 'fr', fm_fr, ...
                                                            'to', fm_to, ...
                                                            'ub', fm_ub );
            end

            pretty_print@mp.dme_branch(obj, dm, section, out_e, mpopt, fd, pp_args);
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = true;
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            fm = pp_args.branch.flow_magnitude;
            rows = find( obj.tab.status & fm.ub ~= 0 & ( ...
                        fm.fr > fm.ub - obj.ctol | ...
                        fm.to > fm.ub - obj.ctol | ...
                        obj.tab.mu_flow_fr_ub > obj.ptol | ...
                        obj.tab.mu_flow_to_ub > obj.ptol ));
        end

        function str = pp_get_title_lim(obj, mpopt, pp_args)
            %
            switch upper(mpopt.opf.flow_lim(1))
                case {'P', '2'}     %% active power
                    str = 'P in MW';
                case 'I'            %% current
                    str = 'I in kA*basekV';
                otherwise           %% apparent power
                    str = 'S in MVA';
            end
            str = sprintf('%s Constraints            (%s)', obj.label, str);
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            switch upper(mpopt.opf.flow_lim(1))
                case {'P', '2'}     %% active power
                    h2 = '   ID      Bus ID    mu_pl_fr   pl_fr    pl_ub    pl_to    mu_pl_to  Bus ID';
                case 'I'            %% current
                    h2 = '   ID      Bus ID    mu_cm_fr   cm_fr    cm_ub    cm_to    mu_cm_to  Bus ID';
                otherwise           %% apparent power
                    h2 = '   ID      Bus ID    mu_sm_fr   sm_fr    sm_ub    sm_to    mu_sm_to  Bus ID';
            end
            h = [ pp_get_headers_lim@mp.dme_shared_opf(obj, dm, out_e, mpopt, pp_args) ...
                {   ' Branch     From        "From" End       Limit       "To" End         To', ...
                    h2, ...
                    '--------  --------  ---------  -------  -------  -------  ---------  --------' } ];
            %%       1234567 123456789 123456.890 12345.78 12345.78 12345.78 123456.890 123456789
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            fm = pp_args.branch.flow_magnitude;

            if fm.ub(k) ~= 0 & (fm.fr(k) > fm.ub(k) - obj.ctol || ...
                    obj.tab.mu_flow_fr_ub(k) > obj.ptol)
                mu_fr = sprintf('%10.3f', obj.tab.mu_flow_fr_ub(k));
            else
                mu_fr = '      -   ';
            end
            if fm.ub(k) ~= 0 & (fm.to(k) > fm.ub(k) - obj.ctol || ...
                    obj.tab.mu_flow_to_ub(k) > obj.ptol)
                mu_to = sprintf('%10.3f', obj.tab.mu_flow_to_ub(k));
            else
                mu_to = '      -   ';
            end

            str = sprintf('%7d %9d %10s %8.2f %8.2f %8.2f %10s %9d', ...
                obj.tab.uid(k), obj.tab.bus_fr(k), mu_fr, fm.fr(k), ...
                fm.ub(k), fm.to(k), mu_to, obj.tab.bus_to(k));
        end
    end     %% methods
end         %% classdef
