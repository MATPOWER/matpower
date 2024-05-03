classdef mp_table_subclass
% mp_table_subclass - Class that acts like a table but isn't one.
%
% Addresses two issues with inheriting from **table** classes (:class:`table`)
% or mp_table).
%
%   1. In MATLAB, :class:`table` is a sealed class, so you cannot inherit
%      from it. You can, however, use a subclass of mp_table, but that can
%      result in the next issue under Octave.
%   2. While nesting of tables works just fine in general, when using mp_table
%      in Octave (at least up through 8.4.0), you cannot nest a subclass of
%      mp_table inside another mp_table object because of this bug:
%      https://savannah.gnu.org/bugs/index.php?65037.
%
% To work around these issues, your "table subclass" can inherit from **this**
% class. An object of this class **isn't** a :class:`table` or mp_table object,
% but rather it **contains** one and attempts to act like one. That is, it
% delegates method calls (currently only those available in mp_table, listed
% below) to the contained table object.
%
% The class of the contained table object is either :class:`table` or mp_table
% and is determined by mp_table_class.
%
% .. admonition:: Limitations
%
%   1. The Octave bug mentioned above also affects tables that inherit from
%      mp_table_subclass. That is, such tables can be nested inside tables
%      of type :class:`table` or mp_table, but not inside tables that are
%      or inherit from mp_table_subclass.
%   2. In MATLAB, when nesting an mp_table_subclass object within another
%      mp_table_subclass object, one cannot use multi-level indexing directly.
%      E.g. If ``T2`` is a variable in ``T1`` and ``x`` is a variable in
%      ``T2``, attempting ``x = T1.T2.x`` will result in an error. The
%      indexing must be done in multiple steps ``T2 = T1.T2; x = T2.x``.
%      Note: This only applies to MATLAB, where the contained table is a
%      :class:`table`. It works just fine in Octave, where the contained
%      table is an :class:`mp_table`.
%
% .. important::
%
%   Since the dot syntax ``T.<var_name>`` is used to access table variables,
%   you must use a functional syntax ``<method>(T,...)``, as opposed to
%   the object-oriented ``T.<method>(...)``, to call methods of this class
%   or subclasses, as with mp_table.
%
% mp.mp_table_subclass Properties:
%   * tab - *(table or mp_table)* contained table object this class emulates
%
% mp.cost_table Methods:
%   * mp_table_subclass - construct object
%   * get_table - return the table stored in :attr:`tab`
%   * set_table - assign a table to :attr:`tab`
%   * istable - true for mp_table objects
%   * size - dimensions of table
%   * isempty - true if table has no columns or no rows
%   * end - used to index last row or variable/column
%   * subsref - indexing a table to retrieve data
%   * subsasgn - indexing a table to assign data
%   * horzcat - concatenate tables horizontally
%   * vertcat - concatenate tables vertically
%   * display - display table contents
%
% See also mp_table, mp_table_class.

%   MATPOWER
%   Copyright (c) 2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tab
    end     %% properties

    %% delegate all mp_table methods to obj.tab
    methods
        function obj = mp_table_subclass(varargin)

            args = varargin;
            if nargin
                %% extract named arguments
                [var_names, row_names, dim_names, args] = ...
                    mp_table.extract_named_args(args);

                %% set default variable names
                nv = length(args);          %% number of variables
                if length(var_names) < nv
                    for k = nv:-1:length(var_names)+1
                        var_names{k} = inputname(k);
                        if isempty(var_names{k})
                            var_names{k} = sprintf('Var%d', k);
                        end
                    end
                end
                args(end+1:end+2) = {'VariableNames', var_names};
                if ~isempty(row_names)
                    args(end+1:end+2) = {'RowNames', row_names};
                end
                if ~isempty(dim_names)
                    args(end+1:end+2) = {'DimensionNames', dim_names};
                end
            end

            table_class = mp_table_class();
            obj.tab = table_class(args{:});
        end

        function tab = get_table(obj)
            %
            % ::
            %
            %   T = get_table(obj)
            tab = obj.tab;
        end

        function obj = set_table(obj, T)
            %
            % ::
            %
            %   set_table(obj, T)
            obj.tab = T;
        end

        function TorF = istable(obj)
            TorF = istable(obj.tab);
        end

        function varargout = size(obj, varargin)
            [varargout{1:nargout}] = size(obj.tab, varargin{:});
        end

        function TorF = isempty(obj)
            TorF = isempty(obj.tab);
        end

        function N = end(obj, k, n)
            N = size(obj.tab, k);
        end

        function n = numArgumentsFromSubscript(obj, varargin)
            n = 1;
        end

        % This is needed to avoid an error when doing T{r1:rN, c} = b;
        % See https://github.com/apjanke/octave-tablicious/issues/80#issuecomment-855198281.
        function n = numel(obj, varargin)
            n = 1;
        end

        function b = subsref(obj, varargin)
            b = subsref(obj.tab, varargin{:});
            if varargin{1}(1).type(1) == '(' && ...
                    (isa(b, 'table') || isa(b, 'mp_table'))
                o = feval(class(obj));
                b = set_table(o, b);
            end
        end

        function obj = subsasgn(obj, varargin)
            for k = 2:length(varargin)
                if isa(varargin{k}, 'mp_table_subclass')
                    varargin{k} = get_table(varargin{k});
                end
            end
            obj.tab = subsasgn(obj.tab, varargin{:});
        end

        function obj = horzcat(obj, varargin)
            args = varargin;
            for k = 1:length(args)
                if isa(args{k}, 'mp_table_subclass')
                    args{k} = args{k}.tab;
                end
            end
            obj.tab = horzcat(obj.tab, args{:});
        end

        function obj = vertcat(obj, varargin)
            args = varargin;
            for k = 1:length(args)
                if isa(args{k}, 'mp_table_subclass')
                    args{k} = args{k}.tab;
                end
            end
            obj.tab = vertcat(obj.tab, args{:});
        end

        function display(obj)
            obj.tab
        end
    end     %% methods
end         %% classdef
