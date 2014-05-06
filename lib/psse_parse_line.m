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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
    if have_fcn('octave')
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
