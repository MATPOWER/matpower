function [data, warns] = psse_parse(records, sections, verbose, rev)
%PSSE_PARSE  Parses the data from a PSS/E RAW data file.
%   DATA = PSSE_PARSE(RECORDS, SECTIONS)
%   DATA = PSSE_PARSE(RECORDS, SECTIONS, VERBOSE)
%   DATA = PSSE_PARSE(RECORDS, SECTIONS, VERBOSE, REV)
%   [DATA, WARNINGS] = PSSE_PARSE(RECORDS, SECTIONS, ...)
%
%   Parses the data from a PSS/E RAW data file (as read by PSSE_READ)
%   into a struct.
%
%   Inputs:
%       RECORDS : cell array of strings, corresponding to the lines
%                 in the RAW file
%       SECTIONS : struct array with indices marking the beginning
%                  and end of each section, and the name of the
%                  section, fields are:
%           first   : index into RECORDS of first line of section
%           last    : index into RECORDS of last line of section
%           name    : name of the section, as extracted from the
%                     END OF ... DATA comments
%       VERBOSE      :  1 (default) to display progress info, 0 otherwise
%       REV          :  (optional) assume the input file is of this
%                       PSS/E revision number, attempts to determine
%                       REV from the file by default
%
%   Output(s):
%       DATA :  a struct with the following fields, each with two
%               sub-fields, 'num' and 'txt' containing the numeric and
%               text data read from the file for the corresponding section
%           id
%           bus
%           load
%           gen
%           shunt
%           branch
%           trans2
%           trans3
%           area
%           twodc
%           swshunt
%       WARNINGS :  cell array of strings containing accumulated
%                   warning messages
%
%   See also PSSE2MPC, PSSE_READ, PSSE_PARSE_SECTION, PSSE_PARSE_LINE

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   Based on mpreadraw.m, written by: Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default args
if nargin < 4
    rev = 0;
    if nargin < 3
        verbose = 0;
    end
end
defaultrev = 23;

%% inititialize section counter
s = 1;
warns = {};

%%-----  case identification data  -----
%% In some files, SBASE is followed by a comment with the REV number such as:
%%  / PSS/E-29.0
%%  / PSS(tm)E-30 RAW
if verbose
    if rev
        fprintf('Forcing interpretation as PSS/E revision %d\n', rev);
    else
        fprintf('Attempting to determine PSS/E revision from content.\n');
    end
    fprintf('Parsing case identification data ...');
end
if rev
    warns{end+1} = sprintf('Conversion explicitly using PSS/E revision %d', rev);
end
[d, c] = psse_parse_line(records{1}, 'dfdfff');
nn = length(d);
data.id.IC = d{1};
if isempty(d{2}) || d{2} <= 0
    error('ERROR: Probable corrupt file, unable to read a valid SBASE value from the first line.');
else
    data.id.SBASE = d{2};
end
if ~isempty(d{3})
    data.id.REV = d{3};
else    %% attempt to extract revision from comment
    tmp = regexp(c, 'PSS(/|\(tm\))E-(?<rev>\d+)', 'tokens');
    if ~isempty(tmp) && size(tmp{1}, 2) == 2
        data.id.REV = str2num(tmp{1}{2});
    else
        data.id.REV = 0;
    end
end
if ~isempty(d{4})
    data.id.XFRRAT = d{4};
else
    data.id.XFRRAT = 0;
end
if ~isempty(d{5})
    data.id.NXFRAT = d{5};
else
    data.id.NXFRAT = 0;
end
if ~isempty(d{6})
    data.id.BASFRQ = d{6};
else
    data.id.BASFRQ = 0;
end
data.id.comment0 = c;
data.id.comment1 = records{2};
data.id.comment2 = records{3};
if verbose
    if rev
        if data.id.REV
            fprintf('.. override detected rev %2d w/%2d ... done.\n', data.id.REV, rev);
        else
            fprintf('...... unknown rev, using rev %2d ... done.\n', rev);
        end
    else
        if data.id.REV
            fprintf('......... rev %2d format detected ... done.\n', data.id.REV);
        else
            fprintf('...... unknown rev, using rev %2d ... done.\n', defaultrev);
        end
    end
end
if ~rev
    if data.id.REV
        rev = data.id.REV;      %% use detected value
    else
        rev = defaultrev;       %% none detected, use default value
        data.id.REV = defaultrev;
        warns{end+1} = sprintf('Unknown REV, using REV %2d format.', defaultrev);
    end
else
    data.id.REV = rev;          %% use override value
end
if isempty(data.id.IC) || data.id.IC ~= 0
    warns{end+1} = sprintf('IC = %d indicates that this may be a change case, rather than base case\n         PSSE2MPC is NOT designed to handle change cases.', data.id.IC);
    if verbose
        fprintf('WARNING : %s\n', warns{end});
    end
end
s = s + 1;

%%-----  bus data  -----
if rev < 24         %% includes load data 
    [data.bus, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'bus', 'ddffffdffsfd');
elseif rev < 31     %% includes fixed shunt data, load separate
    [data.bus, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'bus', 'dsfdffddffd');
else                %% fixed shunt and load data separate
    [data.bus, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'bus', 'dsfddddffff..');
%       'bus', 'dsfddddffffff');
end
s = s + 1;

%%-----  load data  -----
if rev >= 24
    [data.load, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'load', 'd.d..ffffff...');
%       'load', 'dsdddffffffddd');
    s = s + 1;
end

%%-----  fixed shunt data  -----
if rev > 30     %% fixed shunt data is included in bus data for rev <= 30
    [data.shunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'fixed shunt', 'd.dff');
%       'fixed shunt', 'dsdff');
    s = s + 1;
end

%%-----  generator data  -----
[data.gen, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
    'generator', 'd.fffff.f.....d.ff...........');
%   'generator', 'dsfffffdffffffdfffdfdfdfdfsdf');
s = s + 1;

%%-----  branch data  -----
if rev <= 27   %% includes transformer ratio, angle
    [data.branch, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'branch', 'dd.ffffffffffffd');
%       'branch', 'dddffffffffffffd');
elseif rev < 31
    [data.branch, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'branch', 'dd.ffffffffffd');
%       'branch', 'ddsffffffffffdfdfdfdfdf');
else
    [data.branch, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
        'branch', 'dd.ffffffffffd');
%       'branch', 'ddsffffffffffddfdfdfdfdf');
end
s = s + 1;

%%-----  skip transformer adjustment data  -----
if rev <= 27
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'transformer adjustment');
end

%%-----  transformer data  -----
if rev > 27
    %% PSS/E stores two winding and three winding transformer data in the same
    %% section in RAW file. We read in 2 passes, first pass determines the type of
    %% each, second pass reads the data.
    label = 'transformer';
    if ~isempty(sections(s).name) && ~strcmpi(label, sections(s).name)
        if verbose > 1
            fprintf('-----  WARNING:  Expected section labeled: ''%s''\n', upper(label));
            fprintf('-----            Found section labeled:    ''%s''\n', sections(s).name);
        end
    end

    %% Step 1 : Count and collect transformer types
    if verbose
        fprintf('Analyzing transformer types ...');
    end

    %% estimate max number of transformers by number of lines in section
    nt2 = round((sections(s).last - sections(s).first + 1) / 4);
    nt3 = round((sections(s).last - sections(s).first + 1) / 5);

    %% initialize indices for transformer types
    idx2 = zeros(nt2, 1);
    idx3 = zeros(nt3, 1);

    %% set up counters
    i = sections(s).first;      %% initialize record index
    i2 = 0;
    i3 = 0;

    while i <= sections(s).last
        %% determine transformer type
        pat = '[^''",\s/]+\s*(,|\s)\s*[^''",\s/]+\s*(,|\s)\s*([^''",\s/]+)';
        m = regexp(records{i}, pat, 'tokens', 'once');
        if length(m) ~= 3
            disp(m);
            error('m should be length 3');
        end
        if length(m{3}) == 1 && m{3}(1) == '0'  %% two-winding
            i2 = i2 + 1;
            idx2(i2) = i;
            i = i + 4;
        else                                    %% three-winding
            i3 = i3 + 1;
            idx3(i3) = i;
            i = i + 5;
        end
    end
    nt2 = i2;
    nt3 = i3;

    if verbose
        str = sprintf(' %d(%d) two(three)-winding.', nt2, nt3);
        spacers = repmat('.', 1, 36-length(str));
        fprintf('%s %s ... done.\n', spacers, str);
    end

    %% trim index vectors down to size
    idx2 = idx2(1:nt2);
    idx3 = idx3(1:nt3);

    %% parse record 1 (cols 1-20)
    [t2_1, warns] = psse_parse_section(warns, records(idx2), verbose, ...
        '2-winding transformers (1)', 'dd..ddd....d........');
%       '2-winding transformers (1)', 'dddsdddffdsddfdfdfdf');
    [t3_1, warns] = psse_parse_section(warns, records(idx3), verbose, ...
        '3-winding transformers (1)', 'ddd.ddd....d........');
%       '3-winding transformers (1)', 'dddsdddffdsddfdfdfdf');

    %% two-winding
    %% parse record 2 (cols 21-23)
    [t2_2, warns] = psse_parse_section(warns, records(idx2+1), verbose, ...
        '2-winding transformers (2)', 'fff');

    %% parse record 3 (cols 24-39)
    %% parse up to CX1, should warn if CNXA1 is present and non-zero
    [t2_3, warns] = psse_parse_section(warns, records(idx2+2), verbose, ...
        '2-winding transformers (3)', 'ffffff..........');
%       '2-winding transformers (3)', 'ffffffddffffddff');

    %% parse record 4 (cols 40-41)
    [t2_4, warns] = psse_parse_section(warns, records(idx2+3), verbose, ...
        '2-winding transformers (4)', 'ff');

    %% three-winding
    %% parse record 2 (cols 21-31)
    [t3_2, warns] = psse_parse_section(warns, records(idx3+1), verbose, ...
        '3-winding transformers (2)', 'fffffffffff');
%       '3-winding transformers (2)', 'fffffffffff');

    %% parse record 3 (cols 32-47)
    %% parse up to CX1, should warn if CNXA1 is present and non-zero
    [t3_3, warns] = psse_parse_section(warns, records(idx3+2), verbose, ...
        '3-winding transformers (3)', 'ffffff..........');
%       '3-winding transformers (3)', 'ffffffddffffddff');

    %% parse record 4 (cols 48-63)
    %% parse up to CX2
    [t3_4, warns] = psse_parse_section(warns, records(idx3+3), verbose, ...
        '3-winding transformers (4)', 'ffffff..........');
%       '3-winding transformers (4)', 'ffffffddffffddff');

    %% parse record 5 (cols 64-79)
    %% parse up to CX3
    [t3_5, warns] = psse_parse_section(warns, records(idx3+4), verbose, ...
        '3-winding transformers (5)', 'ffffff..........');
%       '3-winding transformers (5)', 'ffffffddffffddff');

    %% assemble two-winding transformer records
    data.trans2.num = [t2_1.num(:, 1:20) t2_2.num(:, 1:3) t2_3.num(:, 1:16) t2_4.num(:, 1:2)];
    data.trans2.txt = [t2_1.txt(:, 1:20) t2_2.txt(:, 1:3) t2_3.txt(:, 1:16) t2_4.txt(:, 1:2)];

    %% assemble three-winding transformer records
    data.trans3.num = [t3_1.num(:, 1:20) t3_2.num(:, 1:11) t3_3.num(:, 1:16) t3_4.num(:, 1:16) t3_5.num(:, 1:16)];
    data.trans3.txt = [t3_1.txt(:, 1:20) t3_2.txt(:, 1:11) t3_3.txt(:, 1:16) t3_4.txt(:, 1:16) t3_5.txt(:, 1:16)];

    % if verbose
    %     fprintf('%s\n', upper(label));
    %     fprintf('%s\n', sections(s).name);
    % end
    s = s + 1;
end

%%-----  area interchange data  -----
[data.area, warns] = psse_parse_section(warns, records, sections, s, verbose, 'area', 'ddffs');
s = s + 1;

%%-----  two-terminal DC transmission line data  -----
label = 'two-terminal DC';
if ~isempty(sections(s).name) && ~strcmpi(label, sections(s).name)
    if verbose > 1
        fprintf('-----  WARNING:  Expected section labeled: ''%s''\n', upper(label));
        fprintf('-----            Found section labeled:    ''%s''\n', sections(s).name);
    end
end
idx = sections(s).first:3:sections(s).last;
if rev < 31
    [dc1, warns] = psse_parse_section(warns, records(idx), verbose, ...
        'two-terminal DC (1)', '.d.ff.......');
%       'two-terminal DC (1)', 'ddffffffsfdf');
    [dc2, warns] = psse_parse_section(warns, records(idx+1), verbose, ...
        'two-terminal DC (2)', 'd.ff.............');
%       'two-terminal DC (2)', 'ddffffffffffdddsf');
    [dc3, warns] = psse_parse_section(warns, records(idx+2), verbose, ...
        'two-terminal DC (3)', 'd.ff.............');
%       'two-terminal DC (3)', 'ddffffffffffdddsf');
else
    [dc1, warns] = psse_parse_section(warns, records(idx), verbose, ...
        'two-terminal DC (1)', '.d.ff.......');
%       'two-terminal DC (1)', 'sdffffffsfdf');
    [dc2, warns] = psse_parse_section(warns, records(idx+1), verbose, ...
        'two-terminal DC (2)', 'd.ff.............');
%       'two-terminal DC (2)', 'ddffffffffffdddDf');
    [dc3, warns] = psse_parse_section(warns, records(idx+2), verbose, ...
        'two-terminal DC (3)', 'd.ff.............');
%       'two-terminal DC (3)', 'ddffffffffffdddDf');
end

%% assemble two-terminal DC transmission line
data.twodc.num = [dc1.num dc2.num dc3.num];
data.twodc.txt = [dc1.txt dc2.txt dc3.txt];
% if verbose
%     fprintf('%s\n', upper(label));
%     fprintf('%s\n', sections(s).name);
% end
s = s + 1;

%%-----  skip voltage source converter data  -----
if rev > 28
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'voltage source converter');
end

%%-----  switched shunt data  -----
if rev < 31
    %% parse up to B1
    if rev <= 27
        [data.swshunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
            'switched shunt', 'd....f');
%           'switched shunt', 'ddffdfdfdfdfdfdfdfdfdf');
    elseif rev <= 29
        [data.swshunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
            'switched shunt', 'd.....f');
%           'switched shunt', 'ddffdsfdfdfdfdfdfdfdfdf');
    else    %%  rev == 30
        [data.swshunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
            'switched shunt', 'd......f');
%           'switched shunt', 'ddffdfsfdfdfdfdfdfdfdfdf');
    end
    s = s + 1;
end

%%-----  skip impedance correction data  -----
[s, warns] = psse_skip_section(warns, sections, s, verbose, 'impedance correction');

%%-----  skip multi-terminal DC data  -----
[s, warns] = psse_skip_section(warns, sections, s, verbose, 'multi-terminal DC');

%%-----  skip multi-section line data  -----
[s, warns] = psse_skip_section(warns, sections, s, verbose, 'multi-section line');

%%-----  skip zone data  -----
[s, warns] = psse_skip_section(warns, sections, s, verbose, 'zone');

%%-----  skip inter-area transfer data  -----
[s, warns] = psse_skip_section(warns, sections, s, verbose, 'inter-area transfer');

%%-----  skip owner data  -----
if rev > 24
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'owner');
end

%%-----  skip FACTS control device data  -----
if rev > 25
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'FACTS control device');
end

%%-----  switched shunt data  -----
if rev > 30
    %% parse up to B1
    if rev < 32
        [data.swshunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
            'switched shunt', 'd......f');
%           'switched shunt', 'ddffdfsfdfdfdfdfdfdfdfdf');
    else
        [data.swshunt, warns] = psse_parse_section(warns, records, sections, s, verbose, ...
            'switched shunt', 'd........f');
%           'switched shunt', 'ddddffdfsfdfdfdfdfdfdfdfdf');
    end
    s = s + 1;
end

%%-----  skip GNE device data  -----
if rev > 31
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'GNE device');
end

%%-----  skip induction machine data  -----
if rev > 32
    [s, warns] = psse_skip_section(warns, sections, s, verbose, 'induction machine');
end

%%-----  check for extra sections  -----
if s <= length(sections)
    warns{end+1} = sprintf('Found %d additional section(s)', length(sections)-s+1);
    if verbose > 1
        fprintf('-----  WARNING:   Found %d additional section(s):\n', length(sections)-s+1);
    end
end
while s <= length(sections)
    n = sections(s).last - sections(s).first + 1;
    if n
        str = sprintf('with %d line(s)', n);
    else
        str = sprintf('(empty)');
    end
    if isempty(sections(s).name)
        warns{end+1} = sprintf('  unlabeled section %s', str);
        if verbose > 1
            fprintf('-----            unlabeled section %s\n', str);
        end
    else
        warns{end+1} = sprintf('  ''%s DATA'' %s', sections(s).name, str);
        if verbose > 1
            fprintf('-----            ''%s DATA'' %s\n', sections(s).name, str);
        end
    end
    s = s + 1;
end



%%---------------------------------------------------------------------
function [s, warns] = psse_skip_section(warns, sections, s, verbose, label)
%PSSE_SKIP_SECTION  Skips over a section without extracting any data.
%   [SIDX, WARNINGS] = PSSE_SKIP_SECTION(WARNINGS, SECTIONS, SIDX, VERBOSE, LABEL)

if s > length(sections)
    if verbose
        spacers = repmat('.', 1, 58-length(label));
        fprintf('No %s data read %s done.\n', label, spacers);
    end
else
    nr = sections(s).last - sections(s).first + 1;
    if nr > 1
        ss = 'lines';
    else
        ss = 'line';
    end
    if nr
        warns{end+1} = sprintf('Skipped %d %s of %s data.', nr, ss, label);
    end
    if ~isempty(sections(s).name) && ~strcmp(upper(label), sections(s).name)
        warns{end+1} = sprintf('Section label mismatch, found ''%s'', expected ''%s''', ...
            sections(s).name, upper(label));
        if verbose
            fprintf('-----  WARNING:  Found section labeled:    ''%s''\n', sections(s).name);
            fprintf('-----            Expected section labeled: ''%s''\n', upper(label));
        end
    end
    if verbose && nr
        spacers = repmat('.', 1, 47-length(ss)-length(label));
        fprintf('Skipping%6d %s of %s data %s done.\n', nr, ss, label, spacers);
    end
    s = s + 1;
end
