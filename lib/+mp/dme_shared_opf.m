classdef (Abstract) dme_shared_opf < handle
% mp.dme_shared_opf - Mixin class for OPF **data model element** objects.
%
% For all elements of mp.data_model_opf, adds shared functionality for
% pretty-printing of **lim** sections.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        ctol    % constraint violation tolerance
        ptol    % shadow price tolerance
    end     %% properties

    methods
        function obj = pp_set_tols_lim(obj, mpopt)
            %
            obj.ctol = mpopt.opf.violation;
            obj.ptol = 1e-4;
        end

        function TorF = pp_have_section_other(obj, section, mpopt, pp_args)
            %
            switch section
                case 'lim'
                    TorF = obj.pp_have_section_lim(mpopt, pp_args);
                otherwise
                    error('mp.dme_shared_opf:pp_have_section_other: unknown section ''%s''', section);
            end
        end

        function rows = pp_rows_other(obj, dm, section, out_e, mpopt, pp_args)
            %
            switch section
                case 'lim'
                    if isempty(obj.ptol)
                        obj.pp_set_tols_lim(mpopt);
                    end
                    rows = obj.pp_rows_lim(dm, out_e, mpopt, pp_args);
                otherwise
                    error('mp.dme_shared_opf:pp_rows_other: unknown section ''%s''', section);
            end
        end

        function h = pp_get_headers_other(obj, dm, section, out_e, mpopt, pp_args)
            %
            switch section
                case 'lim'
                    h = obj.pp_get_headers_lim(dm, out_e, mpopt, pp_args);
                otherwise
                    error('mp.dme_shared_opf:pp_get_headers_other: unknown section ''%s''', section);
            end
        end

        function f = pp_get_footers_other(obj, dm, section, out_e, mpopt, pp_args)
            %
            switch section
                case 'lim'
                    f = obj.pp_get_footers_lim(dm, out_e, mpopt, pp_args);
                otherwise
                    error('mp.dme_shared_opf:pp_get_footers_other: unknown section ''%s''', section);
            end
        end

        function obj = pp_data_other(obj, dm, section, rows, out_e, mpopt, fd, pp_args)
            %
            switch section
                case 'lim'
                    obj.pp_data_lim(dm, rows, out_e, mpopt, fd, pp_args);
                otherwise
                    error('mp.dme_shared_opf:pp_data_other: unknown section ''%s''', section);
            end
        end

        function TorF = pp_have_section_lim(obj, mpopt, pp_args)
            %
            TorF = false;
        end

        function rows = pp_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            if out_e == 2       %% all rows
                rows = -1;
            elseif out_e == 1   %% binding rows
                rows = obj.pp_binding_rows_lim(dm, out_e, mpopt, pp_args);
            else                %% no rows
                rows = 0;
            end
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, pp_args)
            %
            rows = 0;           %% no rows
        end

        function str = pp_get_title_lim(obj, mpopt, pp_args)
            %
            str = sprintf('%s Constraints', obj.label);
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            str = obj.pp_get_title_lim(mpopt, pp_args);
            if isempty(str)
                h = {};
            else
                h = dm.pp_section_label(str);
            end
        end

        function f = pp_get_footers_lim(obj, dm, out_e, mpopt, pp_args)
            %
            f = {};
        end

        function obj = pp_data_lim(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %
            if ~isempty(rows) && rows(1) == -1  %% all rows
                for k = 1:obj.nr
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_lim(dm, k, out_e, mpopt, fd, pp_args));
                end
            else
                for k = 1:length(rows)
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_lim(dm, rows(k), out_e, mpopt, fd, pp_args));
                end
            end
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, pp_args)
            %
            str = '';
        end
    end     %% methods
end         %% classdef
