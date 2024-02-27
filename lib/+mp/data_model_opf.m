classdef data_model_opf < mp.data_model
% mp.data_model_opf - |MATPOWER| **data model** for OPF tasks.
%
% The purpose of this class is to include OPF-specific subclasses for its
% elements and to handle pretty-printing output for **lim** sections.
%
% mp.data_model_opf Methods:
%   * data_model_opf - constructor, assign default data model element classes
%   * pp_flags - add flags for **lim** sections
%   * pp_section_list - append ``'lim'`` tag for **lim** sections to default list
%   * pp_get_headers_other - construct headers for **lim** section headers
%
% See also mp.data_model.

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
        %% constructor
        function obj = data_model_opf()
            % Constructor, assign default data model element classes.
            %
            % Create an empty data model object and assign the default
            % data model element classes, each specific to OPF.
            % ::
            %
            %   dm = mp.data_model_opf()

            %% call parent constructor
            obj@mp.data_model();
            obj.element_classes = ...
                { @mp.dme_bus_opf, @mp.dme_gen_opf, @mp.dme_load_opf, ...
                    @mp.dme_branch_opf, @mp.dme_shunt_opf };
        end

        function [out, add] = pp_flags(obj, mpopt)
            % Add flags for **lim** sections.
            %
            % See :meth:`mp.data_model.pp_flags`.

            %% call parent
            [out, add_] = pp_flags@mp.data_model(obj, mpopt);
            suppress = add_.suppress;

            %% add limit options
            out_s = struct( ...
                'all', 0, ...
                'any', 0, ...
                'elm', add_.s1 );

            if out.all == 1
                out_s.all = 2;
            elseif out.all == -1 && ~suppress
                out_s.all = mpopt.out.lim.all;
            end

            if out_s.all == -1
                out_s.elm.bus     = ~suppress && mpopt.out.lim.v;
                out_s.elm.branch  = ~suppress && mpopt.out.lim.line;
                out_s.elm.gen     = ~suppress && ...
                                    (mpopt.out.lim.pg || mpopt.out.lim.qg);
            else
                out_s.elm.bus     = out_s.all;
                out_s.elm.branch  = out_s.all;
                out_s.elm.gen     = out_s.all;
            end
            if isfield(mpopt.out.lim, 'elm')
                out_s.elm = nested_struct_copy(out_s.elm, mpopt.out.lim.elm);
            end
            out_s.any = any(cell2mat( ...
                    cellfun(@double, struct2cell(out_s.elm), ...
                            'UniformOutput', false) ));

            out.sec.lim = out_s;
            out.any = out.any || out_s.any;

            %% return additional data
            if nargout > 1
                add = add_;
            end
        end

        function sections = pp_section_list(obj, out)
            % Append ``'lim'`` tag for **lim** section  to default list.
            %
            % See :meth:`mp.data_model.pp_section_list`.

            sections = pp_section_list@mp.data_model(obj, out);
            sections{end+1} = 'lim';
        end

        function h = pp_get_headers_other(obj, section, out_s, mpopt)
            % Construct pretty printed lines for **lim** section headers.
            %
            % See :meth:`mp.data_model.pp_get_headers_other`.

            switch section
                case 'lim'
                    h = {};
                otherwise
                    error('mp.data_model_opf:pp_get_headers_other: unknown section ''%s''', section);
            end
        end
    end     %% methods
end         %% classdef
