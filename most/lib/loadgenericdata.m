function var = loadgenericdata(varfile, vartype, fields, varname, args)
%LOADGENERICDATA Loads data of a specified type from file or data structure.
%
%   VAR = LOADGENERICDATA(VARFILE, VARTYPE)
%   VAR = LOADGENERICDATA(VARFILE, VARTYPE, FIELDS)
%   VAR = LOADGENERICDATA(VARFILE, VARTYPE, FIELDS, VARNAME)
%   VAR = LOADGENERICDATA(VARFILE, VARTYPE, FIELDS, VARNAME, ARGS)
%
%   Loads data from a variable or M-file or MAT-file and checks that it
%   matches a specified type.
%
%   Inputs:
%       VARFILE : Variable containing the data structure or a string
%                 containing the name of a function M-file or a MAT-file
%                 on the MATLAB path. If no file extension is provided
%                 it will attempt to load a MAT-file with that name and,
%                 if not found, will call a function by that name to
%                 get the data. The function M-file should return a
%                 single argument containing the data. A MAT-file should
%                 either contain a single variable with the desired data
%                 or provide the variable name in VARNAME.
%       VARTYPE : String or cell array of strings with, in order of
%                 priority, the data structure type to be returned.
%                 Valid values are 'struct', 'cell' and 'array'.
%       FIELDS :  (optional) String or cell array of strings containing
%                 a list of required fields in case the VARTYPE is struct.
%                 If a required field is missing it will throw an error.
%       VARNAME   (optional) String containing the name of the variable to
%                 extract when loading a MAT-file. If not provided, the
%                 default is to extract the 1st variable, regardless of name.
%       ARGS :    (optional) Scalar or cell array of values that are passed
%                 as input function arguments to VARFILE if it is an M-file.
%
%   Output:
%       VAR : Returned data structure of the first matching VARTYPE.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Daniel Munoz-Alvarez and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if isempty(varfile)
    var = [];
    return;
end
if nargin < 5
    args = [];
    if nargin < 4
        varname = '';
        if nargin < 3
            fields = {};
        end
    end
end

if ischar(varfile)      %% it's a file name
    %% check for explicit extension
    l = length(varfile);
    if l > 2
        if strcmp(varfile(l-1:l), '.m')
            rootname = varfile(1:l-2);
            extension = '.m';
        elseif l > 4
            if strcmp(varfile(l-3:l), '.mat')
                rootname = varfile(1:l-4);
                extension = '.mat';
            end
        end
    end

    %% set extension if not specified explicitly
    if ~exist('rootname', 'var')
        rootname = varfile;
        if exist([rootname '.mat'], 'file') == 2
            extension = '.mat';
        elseif exist([rootname '.m'], 'file') == 2
            extension = '.m';
        else
            error('loadgenericdata: No file named ''%s.mat'' or ''%s.m'' found in MATLAB path.', rootname, rootname);
        end
    end

    %% attempt to read file
    if strcmp(extension,'.mat')         %% from MAT file
        s = load(rootname);
        if isempty(varname)
            f = fieldnames(s);
            var = s.(f{1});             %% use first variable if not specified
        elseif isfield(s, varname)
            var = s.(varname);          %% use specified variable
        else
            error('loadgenericdata: No variable named ''%s'' in ''%s.mat''', varname, rootname);
        end
    elseif strcmp(extension,'.m')       %% from M file
        if isempty(args)
            var = feval(rootname);
        elseif iscell(args)
            var = feval(rootname, args{:});
        else
            var = feval(rootname, args);
        end
    end
else                    %% it's a data structure, not a file name
    var = varfile;
end

%% check the provided types
if ~iscell(vartype)
    vartype = {vartype};    %% convert string to single element cell array
end
if ~iscell(fields) && ~isempty(fields)
    fields = {fields};      %% convert string to single element cell array
end
done = false;
for i = length(vartype)
    switch vartype{i}
        case 'struct'
            if isstruct(var)
                for j = 1:length(fields)
                    if ~ischar(fields{j})
                        error('loadgenericdata: FIELDS input must be a string or cell array of strings.')
                    end
                    if ~isfield(var, fields{j})
                        error('loadgenericdata:missingfield', 'loadgenericdata: field ''%s'' missing from struct data.', fields{j});
                    end
                end
                done = true; break;     %% struct OK
            end

        case 'cell'
            if iscell(var)
                done = true; break;     %% cell array OK
            end

        case 'array'
            if isnumeric(var)
                done = true; break;     %% array OK
            end
        otherwise
            error('loadgenericdata: ''%s'' is not a valid value for VARTYPE', vartype{i});
    end
end
if ~done
    error('loadgenericdata: data not found in any of the specified formats');
end
