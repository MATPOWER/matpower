classdef mp_table
%MP_TABLE  Very basic TABLE class for Octave or older Matlab.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties (Access=protected)
        Properties = struct('VariableNames',  {{}}, ...
                            'VariableValues', {{}}, ...
                            'RowNames',       {{}}, ...
                            'DimensionNames', {{}}    );
    end     %% properties

    methods
        function obj = mp_table(varargin)
            args = varargin;

            if nargin
                %% extract named arguments
                [var_names, row_names, dim_names, args] = ...
                    extract_named_args(obj, args);

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

                %% set default dimension names
                if isempty(dim_names)
                    dim_names = {'Row', 'Variables'};
                end

                %% check dimensions
                if length(var_names) > nv
                    warning('mp_table: ignoring %d extra VariableNames provided', length(var_names)-nv);
                end
                nr = size(args{1}, 1);      %% number of rows
                for k = 2:nv
                    if size(args{k}, 1) ~= nr
                        error('mp_table: variable %d expected to have %d rows', k, nr);
                    end
                end

                %% set table properties
                obj.Properties.VariableNames = var_names;
                obj.Properties.DimensionNames = dim_names;
                obj.Properties.RowNames = row_names;
                obj.Properties.VariableValues = args;
                obj.Properties.Map = cell2struct(num2cell(1:nv), var_names, 2);
            end
        end

        function new_obj = copy(obj)
            new_obj = eval(class(obj));  %% create new object

            %% make copies of properties
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end

            %% make copies of values
            for k = 1:length(obj.Properties.VariableValues)
                if isobject(obj.Properties.VariableValues{k}) && ...
                        ismethod(obj.Properties.VariableValues{k}, 'copy')
                    new_obj.Properties.VariableValues{k} = ...
                        copy(obj.Properties.VariableValues{k});
                end
            end
        end

        %% not really useable until Octave implements istable()
        function TorF = istable(obj)
            TorF = true;
        end

        function varargout = size(obj, dim)
            varargout = cell(1, nargout);
            w = length(obj.Properties.VariableValues);
            if w
                h = size(obj.Properties.VariableValues{1}, 1);
            else
                h = 0;
            end
            if nargin == 2
                if dim == 1
                  varargout{1} = h;
                elseif dim == 2
                  varargout{1} = w;
                else
                  varargout{1} = 1;
                end
            elseif nargout == 0 || nargout == 1
                varargout{1} = [h, w];
            else
                varargout{1} = h;
                varargout{2} = w;
                [varargout{3:end}] = deal(1);
            end
        end

        function TorF = isempty(obj)
            TorF = prod(size(obj)) == 0;
        end

        function N = end(obj, k, n)
            N = size(obj, k);
        end

        function n = numArgumentsFromSubscript(obj, s, indexingContext)
            n = 1;
        end

        function n = numel(obj, varargin)
            %% needed to avoid error when doing T{r1:rN, c} = b;
            %% see https://github.com/apjanke/octave-tablicious/issues/80#issuecomment-855198281
            n = 1;
        end

        function b = subsref(obj, s)
            switch s(1).type
                case '.'
                    if strcmp(s(1).subs, 'Properties')
                        b = obj.Properties;
                    else
                        p = obj.Properties;
                        b = p.VariableValues{p.Map.(s(1).subs)};
                    end
                case '()'
                    r = s(1).subs{1};
                    if length(s(1).subs) > 1
                        c = s(1).subs{2};
                    else
                        c = ':';
                    end
                    var_names = obj.Properties.VariableNames(1,c);
                    %% Octave has issues with this ...
                    % var_vals = cellfun(@(x)x(r,:), ...
                    %     obj.Properties.VariableValues(c), 'UniformOutput', 0);
                    %% ... so we use explicit subsref instead ...
                    ss = struct('type', '()', 'subs', {{r, ':'}});
                    var_vals = cellfun(@(x)subsref(x, ss), ...
                        obj.Properties.VariableValues(c), 'UniformOutput', 0);
                    dim_names = obj.Properties.DimensionNames;
                    if isempty(obj.Properties.RowNames)
                        row_names = {};
                    else
                        row_names = obj.Properties.RowNames(r);
                    end
                    b = mp_table(var_vals{:}, ...
                        'VariableNames', var_names, 'RowNames', row_names, ...
                        'DimensionNames', dim_names);
                case '{}'
                    a = horzcat(obj.Properties.VariableValues{s(1).subs{2}});
                    b = a(s(1).subs{1}, :);
            end
            if length(s) > 1    %% recurse
                b = subsref(b, s(2:end));
            end
        end

        function obj = subsasgn(obj, s, b)
            R = length(s) > 1;      %% need to recurse
            switch s(1).type
                case '.'
                    if strcmp(s(1).subs, 'Properties')
                        if R
                            obj.Properties = ...
                                subsasgn(obj.Properties, s(2:end), b);
                        else
                            obj.Properties = b;
                        end
                    else
                        if isfield(obj.Properties.Map, s(1).subs)
                            idx = obj.Properties.Map.(s(1).subs);
                            if R
                                obj.Properties.VariableValues{idx} = ...
                                    subsasgn(obj.Properties.VariableValues{idx}, ...
                                        s(2:end), b);
                            else
                                obj.Properties.VariableValues{idx} = b;
                            end
                        else    %% new variable
                            if R
                                error('mp_table/subsagn: cannot index new variable ''%s''', s(1).subs);
                            end
                            %% check size compatibility
                            if size(b, 1) ~= size(obj, 1)
                                error('mp_table/subsagn: number of rows (%d) in new variable ''%s'' does not match height of table (%d)', size(b, 1), s(1).subs, size(obj, 1));
                            end
                            %% add variable
                            idx = size(obj, 2) + 1;
                            obj.Properties.VariableNames{idx} = s(1).subs;
                            obj.Properties.Map.(s(1).subs) = idx;
                            obj.Properties.VariableValues{idx} = b;
                        end
                    end
                case '()'
                    nv = size(b, 2);
                    idx = s(1).subs{2};
                    if ischar(idx)
                        if length(idx) == 1 && idx == ':'
                            idx = [1:nv];
                        else        %% for 'a:b' syntax as a char array
                            idx = eval(['[' idx ']'])
                        end
                    end
                    if nv
                        ss = struct('type', '{}', 'subs', {{':', 0}});
                        for k = 1:nv
                            ss(1).subs{2} = k;
                            obj.Properties.VariableValues{idx(k)}(s(1).subs{1}, :) = ...
                                subsref(b, ss);     %% = b{:, k};
                        end
                    else    %% delete variable
                        obj.Properties.VariableNames(idx) = [];
                        obj.Properties.VariableValues(idx) = [];
                        %% rebuild map
                        nv = size(obj, 2);
                        obj.Properties.Map = cell2struct(num2cell(1:nv), obj.Properties.VariableNames, 2);
                    end
                    if R
                        error('mp_table/subasgn: sub-indexing following {}-indexing not supported');
                    end
                case '{}'
                    idx = s(1).subs{2};
                    if ischar(idx)
                        if length(idx) == 1 && idx == ':'
                            idx = [1:size(obj, 2)];
                        else        %% for 'a:b' syntax as a char array
                            idx = eval(['[' idx ']'])
                        end
                    end
                    j = 1;  %% col index in b
                    for k = 1:length(idx)
                        %% get num cols for var k
                        nc = size(obj.Properties.VariableValues{idx(k)}, 2);
                        if nc > 1
                            obj.Properties.VariableValues{idx(k)}(s(1).subs{1}, :) = b(:, j:j+nc-1);
                        else
                            obj.Properties.VariableValues{idx(k)}(s(1).subs{1}, :) = b(:, j);
                        end
                        j = j+nc;
                    end
                    if R
                        error('mp_table/subasgn: sub-indexing following ()-indexing not supported');
                    end
            end
        end

        function obj = horzcat(obj, varargin)
            [nr, nv] = size(obj);
            nvn = length(obj.Properties.VariableNames); %% num variable names
            nrn = length(obj.Properties.RowNames);      %% num row names

            %% for each table to append
            for j = 1:length(varargin)
                Tj = varargin{j};       %% table being appended
                [nrj, nvj] = size(Tj);
                nrnj = length(Tj.Properties.RowNames);      %% num row names Tj

                %% check for consistent number of rows
                if nrj ~= nr
                    error('mp_table/horzcat: number of rows in arg %d (%d) does not match that in arg 1 (%d)', ...
                        j+1, nrj, nr);
                end

                %% check for row name mismatch
                if nrn          %% obj has row names
                    if nrnj && ~isequal(obj.Properties.RowNames, ...
                                         Tj.Properties.RowNames)
                        for k = 1:nrn
                            if ~strcmp(obj.Properties.RowNames{k}, ...
                                        Tj.Properties.RowNames{k})
                                error('mp_table/horzcat: name of row %d in arg %d (%s) does not match the one in arg 1 (%s)', ...
                                    k, j+1, Tj.Properties.RowNames{k}, ...
                                        obj.Properties.RowNames{k});
                            end
                        end
                    end
                elseif nrnj     %% Tj has row names, but obj does not (use them)
                    obj.Properties.RowNames = Tj.Properties.RowNames;
                    nrn = nrnj;
                end
                var_names = Tj.Properties.VariableNames;

                %% stack variable names and check for duplicates
                obj.Properties.VariableNames = ...
                    [obj.Properties.VariableNames, var_names];
                nvn = nvn + nvj;    %% update num variable names
                if length(unique(obj.Properties.VariableNames)) < nvn
                    [~, ia, ~] = unique(obj.Properties.VariableNames);
                    k = find(ismember(1:nvn, ia));
                    error('mp_table/horzcat: duplicate variable name ''%s''', ...
                        obj.Properties.VariableNames{k(1)});
                end

                %% horzcat var values
                obj.Properties.VariableValues = ...
                    [ obj.Properties.VariableValues, ...
                       Tj.Properties.VariableValues ];
                nv = length(obj.Properties.VariableValues);

                %% update Properties.Map
                obj.Properties.Map = cell2struct(num2cell(1:nv), ...
                    obj.Properties.VariableNames, 2);
            end
        end

        function obj = vertcat(obj, varargin)
            [nr, nv] = size(obj);
            nvn = length(obj.Properties.VariableNames); %% num variable names
            nrn = length(obj.Properties.RowNames);      %% num row names

            %% for each table to append
            for j = 1:length(varargin)
                Tj = varargin{j};       %% table being appended
                [nrj, nvj] = size(Tj);
                nvnj = length(Tj.Properties.VariableNames); %% num var names Tj

                %% check for consistent number of variables
                if nvj ~= nv
                    error('mp_table/vertcat: number of variables in arg %d (%d) does not match that in arg 1 (%d)', ...
                        j+1, nvj, nv);
                end

                %% check for variable name mismatch
                if ~isequal(obj.Properties.VariableNames, ...
                             Tj.Properties.VariableNames)
                    for k = 1:nvn
                        if ~strcmp(obj.Properties.VariableNames{k}, ...
                                    Tj.Properties.VariableNames{k})
                            error('mp_table/vertcat: name of variable %d in arg %d (%s) does not match the one in arg 1 (%s)', ...
                                k, j+1, Tj.Properties.VariableNames{k}, ...
                                    obj.Properties.VariableNames{k});
                        end
                    end
                end

                %% create missing row names if any exist
                row_names = Tj.Properties.RowNames;
                if nrn                      %% obj has row names, Tj does not
                    if length(row_names) == 0
                        %% create names for rows of Tj
                        row_names = cell(nrj, 1);
                        for kk = 1:nrj
                            row_names{kk} = sprintf('Row%d', kk+nrn);
                        end
                    end
                elseif length(row_names)    %% Tj has row names, obj does not
                    %% create names for rows of obj
                    obj.Properties.RowNames = cell(nr, 1);
                    for kk = 1:nr
                        obj.Properties.RowNames{kk} = sprintf('Row%d', kk);
                    end
                end
                nrn = length(obj.Properties.RowNames);

                %% stack row names and check for duplicates
                if nrn
                    obj.Properties.RowNames = ...
                        [obj.Properties.RowNames; row_names];
                    nrn = nrn + nrj;    %% update num row names
                    if length(unique(obj.Properties.RowNames)) < nrn
                        [~, ia, ~] = unique(obj.Properties.RowNames);
                        k = find(ismember(1:nrn, ia));
                        error('mp_table/horzcat: duplicate row name ''%s''', ...
                            obj.Properties.RowNames{k(1)});
                    end
                end

                %% for each variable
                for k = 1:nv
                    %% vertcat var values
                    obj.Properties.VariableValues{k} = ...
                        [ obj.Properties.VariableValues{k}; ...
                           Tj.Properties.VariableValues{k} ];
                end
            end
        end

        function [var_names, row_names, dim_names, args] = ...
                extract_named_args(obj, args)
            var_names = {};
            row_names = {};
            dim_names = {};
            for k = length(args)-1:-1:1
                if ischar(args{k})
                    switch args{k}
                        case 'VariableNames'
                            var_names = args{k+1};
                            args(k:k+1) = [];
                        case 'RowNames'
                            row_names = args{k+1};
                            args(k:k+1) = [];
                        case 'DimensionNames'
                            dim_names = args{k+1};
                            args(k:k+1) = [];
                    end
                end
            end
        end

        function display(obj)
            allrows = 1;

            [nr, nv] = size(obj);

            %% rows to display
            max_nr = 25;
            first_last = 10;
            if nr <= max_nr || allrows
                rr = [1:nr]';
            else
                rr = [1:first_last nr-first_last+1:nr]';
            end

            %% var names
            vn = obj.Properties.VariableNames;
            rn = obj.Properties.RowNames;

            %% col_widths
            if isempty(rn)
                cwr = 0;
            else
                cwr = max(cellfun(@length, rn));    %% col width for row names
            end
            for v = nv:-1:1
                cw(v) = length(vn{v});  %% actual column width
                cwv(v) = 0;             %% column width for displayed values
                vv = obj.Properties.VariableValues{v};
                if size(vv, 2) > 1 || isa(vv, 'mp_table') || isa(vv, 'table')
                    s = sprintf('1x%d %s' , size(vv, 2), class(vv));
                    len_s = length(s);
                    if len_s > cw(v);
                        cw(v) = len_s;
                    end
                    if len_s > cwv(v);
                        cwv(v) = len_s;
                    end
                else
                    if iscell(vv)
                        for r = 1:length(rr)
                            s = sprintf('''%s''', vv{rr(r)});
                            len_s = length(s);
                            if len_s > cw(v);
                                cw(v) = len_s;
                            end
                            if len_s > cwv(v);
                                cwv(v) = len_s;
                            end
                        end
                    else
                        for r = 1:length(rr)
                            if isempty(vv)
                                s = '[]';
                            else
                                s = sprintf('%g', vv(rr(r), 1));
                            end
                            len_s = length(s);
                            if len_s > cw(v);
                                cw(v) = len_s;
                            end
                            if len_s > cwv(v);
                                cwv(v) = len_s;
                            end
                        end
                    end
                end
                spc{v} = repmat(' ', 1, floor((cw(v)-cwv(v))/2));
            end

            fprintf('\n %dx%d %s\n\n', nr, nv, class(obj));

            %% display var names & lines
            fprintf('  %s', repmat(' ', 1, cwr));
            for v = 1:nv
                fmt = sprintf('  %%-%ds', cw(v));
                name = sprintf('%s', vn{v});
                spcs = repmat(' ', 1, floor((cw(v)-length(name))/2));
                fprintf(fmt, [spcs name]);
            end
            fprintf('\n');
            fprintf('  %s', repmat(' ', 1, cwr));
            for v = 1:nv
                s = repmat('-', 1, cw(v));
                fprintf('  %s', s);
            end
            fprintf('\n');

            %% display var values
            for r = 1:length(rr)
                %% row names
                if isempty(rn)
                    fprintf('  ');
                else
                    fmt = sprintf('  %%-%ds', cwr);
                    fprintf(fmt, rn{rr(r)});
                end

                %% variable values
                for v = 1:nv
                    vv = obj.Properties.VariableValues{v};
                    if size(vv, 2) > 1 || isa(vv, 'mp_table') || isa(vv, 'table')
                        val = sprintf('%s1x%d %s', spc{v}, size(vv, 2), class(vv));
                        fmt = sprintf('  %%-%ds', cw(v));
                    else
                        if iscell(vv)
                            val = sprintf('%s''%s''', spc{v}, vv{rr(r), 1});
                            fmt = sprintf('  %%-%ds', cw(v));
                        else
                            if isempty(vv)
                                val = sprintf('[]%s', spc{v});
                            else
                                val = sprintf('%g%s', vv(rr(r), 1), spc{v});
                            end
                            fmt = sprintf('  %%%ds', cw(v));
                        end
                    end
                    fprintf(fmt, val);
                end
                fprintf('\n');

                %% print continuation dots if not displaying all rows
                if nr > length(rr) && r == first_last
                    for k = 1:3
                        for v = 1:nv
                            spcs = repmat(' ', 1, floor(cwv(v)/2));
                            if iscell(vv)
                                dot = sprintf('%s%s.', spc{v}, spcs);
                                fmt = sprintf('  %%-%ds', cw(v));
                            else
                                dot = sprintf('.%s%s', spcs, spc{v});
                                fmt = sprintf('  %%%ds', cw(v));
                            end
                            fprintf(fmt, dot);
                        end
                        fprintf('\n');
                    end
                end
            end
            fprintf('\n');
        end
    end     %% methods
end         %% classdef
