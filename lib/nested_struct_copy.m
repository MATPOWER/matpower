function d = nested_struct_copy(d, s, opt, parent)
%NESTED_STRUCT_COPY  Copies values from one nested struct to another.
%
%   DS = NESTED_STRUCT_COPY(D, S)
%   DS = NESTED_STRUCT_COPY(D, S, OPT)
%   DS = NESTED_STRUCT_COPY(D, S, OPT, PARENT)
%
%   DS = NESTED_STRUCT_COPY(D, S) copies values from a source struct S to
%   a destination struct D in a nested, recursive manner. That is, the value
%   of each field in S is copied directly to the corresponding field in D,
%   unless that value is itself a struct, in which case the copy is done
%   via a recursive call to NESTED_STRUCT_COPY.
%
%   Inputs:
%       D - the destination struct that values are copied to
%       S - the source struct containing the values to be copied from
%       OPT - (optional) options struct controlling copy behavior, with fields:
%           check - check that field name is valid, by looking for it in
%                   OPT.valid_fields (defaults to D), before copying
%                0 - (default) do not do any field name checking
%                1 - fatal error if S contains an invalid field name
%               -1 - skip any invalid fields in S
%           copy_mode - how to handle assignment of fields that are structs
%               ''          - (default) recursive call to nested_struct_copy()
%               '='         - direct assignment, D.<field> = S.<field>
%               @<function> - pointer to a function to be called with field
%                             from S, returning field to assign to D
%                             D.<field> = <function>(S.<field>)
%           valid_fields - struct containing, the heirarchy of all of (and
%                       only) the valid field names (field values are ignored)
%           exceptions - a struct array, with the following fields, defining
%                       exceptions to the top-level options
%               name        - name (can be multi-level) of field to which
%                             exception applies
%               check       - same as OPT.check,     only for specified field
%               copy_mode   - same as OPT.copy_mode, only for specified field
%               valid_fields- same as OPT.valid_fields,   for specified field
%       PARENT - cell array of parent field names used by NESTED_STRUCT_COPY
%               with recursive calls to allow checking of multi-level field
%               field names in exceptions, e.g. when called recursively to
%               assign the field S.info.address.home the value of PARENT would
%               be {'info', 'address'}.
%
%   Output:
%       DS - the combined destination struct 
%
%   Examples:
%       See T_NESTED_STRUCT_COPY (t_nested_struct_copy.m).
%
%   TO DO:  Finish example.
%           Implement an error that passes back the full field string of
%           an invalid field so that mpoption can refer to it as option foo.

%   MATPOWER
%   Copyright (c) 2013-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

DEBUG = 0;

%% set default input args
if nargin < 4
    parent = {};
    if nargin < 3
        opt = struct;
    end
end

%% set up options
if isfield(opt, 'check')
    check = opt.check;
else
    check = 0;
end
if isfield(opt, 'copy_mode')
    copy_mode = opt.copy_mode;
else
    copy_mode = '';
end
if isfield(opt, 'valid_fields')
    valid_fields = opt.valid_fields;
else
    valid_fields = d;
end
if isfield(opt, 'exceptions')
    exceptions = opt.exceptions;
else
    exceptions = struct('name', {});
end

%% form parent string
if DEBUG, fprintf('nested_struct_copy() : parent = %s\n', strjoin(parent, '.')); end
if nargin > 3 && ~isempty(parent)
    pl = length(parent);
    tmp = cell(2, pl);
    tmp(1,:) = parent;
    tmp(2,1:pl-1) = {'.'};
    tmp(2,pl) = {''};
    parentstr = [tmp{:}];
else
    parentstr = '';
end

%% process fields
fields = fieldnames(s);
for f = 1:length(fields)
    ff = fields{f};
    
    %% form full field name
    if isempty(parentstr)
        str = ff;
    else
        str = [parentstr '.' ff];
    end
    
    %% field doesn't exist in valid_fields
    if ~isfield(valid_fields, ff)
        if check > 0        %% throw an error
            error('nested_struct_copy: ''%s'' is not a valid field name', str);
        elseif check < 0    %% skip to next field
            continue;
        end
    end
    
    ck = check;
    cm = copy_mode;
    vf = valid_fields;
    
    %% look for an exception that matches this field
    if isempty(exceptions)
        k = [];
    else
        k = find(strcmp(str, {exceptions.name}'), 1);
        if ~isempty(k)
            if isfield(exceptions, 'copy_mode') && ...
                    ( ischar(exceptions(k).copy_mode) || ...
                     ~isempty(exceptions(k).copy_mode) )
                cm = exceptions(k).copy_mode;
            end
            if isfield(exceptions, 'check') && ...
                    ~isempty(exceptions(k).check)
                ck = exceptions(k).check;
            end
            if isfield(exceptions, 'valid_fields') && ...
                    ~isempty(exceptions(k).valid_fields)
                vf = exceptions(k).valid_fields;
            end
        end
    end

    %% copy the field
    if strcmp(class(cm), 'function_handle')
        %% assign via function handle
        d.(ff) = cm(s.(ff));
    elseif ~isstruct(s.(ff)) || (ischar(cm) && strcmp(cm, '=')) || ...
            (isfield(d, ff) && ~isstruct(d.(ff)))
        %% non-struct OR struct with cm == '=' OR struct to non-struct copy
        %% assign directly
        d.(ff) = s.(ff);
    elseif isstruct(s.(ff)) && isempty(cm)
        %% assign via recursive call to nested_struct_copy()
        if isfield(vf, ff)
            newvf = vf.(ff);    %% use sub-field of valid_fields if present
        else
            newvf = struct;
        end
        newopt = struct( ...
            'check',        ck, ...
            'copy_mode',    cm, ...
            'valid_fields', newvf, ...
            'exceptions',   exceptions ...
        );
        if ~isfield(d, ff)  %% create field if it doesn't exist in d
            d.(ff) = struct;
        end
        ss = s.(ff);
        if length(ss) > 1
            d.(ff) = ss;
        else
            d.(ff) = nested_struct_copy(d.(ff), ss, newopt, {parent{:}, ff});
        end
    else
        error('nested_struct_copy: OPT.copy_mode must be '''', ''='', or a function handle\n');
    end
end
