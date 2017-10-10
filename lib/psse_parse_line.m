function [data, comment] = psse_parse_line(str, t)
%PSSE_PARSE_LINE  Reads and parses a single line from a PSS/E RAW data file
%   [DATA, COMMENT] = PSSE_PARSE_LINE(FID)
%   [DATA, COMMENT] = PSSE_PARSE_LINE(FID, TEMPLATE)
%   [DATA, COMMENT] = PSSE_PARSE_LINE(STR)
%   [DATA, COMMENT] = PSSE_PARSE_LINE(STR, TEMPLATE)
%
%   Parses a single line from a PSS/E RAW data file, either directly read
%   from the file, or passed as a string.
%
%   Inputs:
%       FID :       (optional) file id of file from which to read the line
%       STR :       string containing the line to be parsed
%       TEMPLATE :  (optional) string of characters indicating how to
%                   interpret the type of the corresponding column, options
%                   are as follows:
%               d, f or g : integer floating point number to be converted
%                   via SSCANF with %d, %f or %g, respectively.
%               D, F or G : integer floating point number, possibly enclosed
%                   in single or double quotes, to be converted via
%                   SSCANF with %d, %f or %g, respectively.
%               c or s : character or string, possibly enclosed in single
%                   or double quotes, which are stripped from the string
%           Note:   Data columns in STR that have no valid corresponding
%                   entry in TEMPLATE (beyond end of TEMPLATE, or a character
%                   other than those listed, e.g. '.') are returned as a
%                   string with no conversion. TEMPLATE entries for which
%                   there is no corresponding column are returned as NaN or
%                   empty string, depending on the type.
%   Outputs:
%       DATA :      a cell array whose elements contain the contents of
%                   the corresponding column in the data, converted
%                   according to the TEMPLATE.
%       COMMENT :   (optional) possible comment at the end of the line

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% read the line
if ischar(str)
    ln = str;
else
    ln = fgets(str);
end

if ischar(ln)
    %% parse the line
    delim = '\s*(,|\s)\s*';     %% general delimiter
    repeatdelim = '\s*,\s*|\t'; %% delimiter that allows repeated delimiters
    non_quote_field = '[^''",\s/]+';
    single_quote_field = '''([^'']|'''')*''';
    double_quote_field = '"([^"]|"")*"';
    any_field = sprintf('(?<col>%s|%s|%s)', non_quote_field, single_quote_field, double_quote_field);
    pat = sprintf('%s%s|%s|%s|(?<comment>/.*)?', any_field, delim, repeatdelim, any_field);
    % pat = '(?<col>[^''",\s/]+|''([^'']|'''')*''|"([^"]|"")*")\s*(,|\s)\s*|\s*,\s*|\t|(?<col>[^''",\s/]+|''([^'']|'''')*''|"([^"]|"")*")|(?<comment>/.*)?';
    n = regexp(ln, pat, 'names');

    %% extract data
    if have_fcn('octave') && have_fcn('octave', 'vnum') < 4.003
        nc = length(n.col);
        if nc && isempty(n.col{nc})
            data = n.col(1:nc-1);
            if ~isempty(n.comment{nc})
                comment = strtrim(n.comment{nc}(2:end));
            else
                comment = '';
            end
            nc = nc - 1;
        else
            data = n.col;
            comment = '';
        end
    else
        nc = length(n);
        if nc && isempty(n(nc).col)
            [data{1:nc-1}] = deal(n(1:nc-1).col);
            if ~isempty(n(nc).comment)
                comment = strtrim(n(nc).comment(2:end));
            else
                comment = '';
            end
            nc = nc - 1;
        else
            [data{1:nc}] = deal(n.col);
            comment = '';
        end
    end

    %% check for section end and Q record
    if length(data) == 1 && length(data{1}) == 1
        if data{1}(1) == '0'
            data{1} = 0;
        elseif data{1}(1) == 'Q'
            data = {};
        end
    end
else
    data = {};
    comment = '';
end

%% clean/convert data (convert numeric, strip quotes from strings)
if nargin > 1 && ~isempty(t) && ~isempty(data) && ...
        (length(data) ~= 1 || ~isnumeric(data{1}) || data{1} ~= 0)
    nt = length(t);
    for k = 1:min(nc, nt)
        switch t(k)
            case {'D', 'F', 'G', 's', 'c'}
                if ~isempty(data{k}) && (data{k}(1) == '''' || data{k}(1) == '"')
                    data{k} = data{k}(2:end-1);
                end
            % otherwise             %% do nothing
        end
        switch upper(t(k))
            case {'D', 'F', 'G'}
                data{k} = sscanf(data{k}, ['%' lower(t(k))]);
            % otherwise             %% do nothing
        end
    end
    if nc < nt
        data(nc+1:nt) = cell(1,nt-nc);
    end
end
