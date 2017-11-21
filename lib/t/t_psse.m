function t_psse(quiet)
%T_PSSE  Tests for PSSE2MPC and related functions.

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

num_tests = 142;

t_begin(num_tests, quiet);

raw = 't_psse_case.raw';
case_n = 't_psse_case%d';
casefile = 't_case9_save2psse';
if quiet
    verbose = 0;
else
    verbose = 0;
end
if have_fcn('octave')
    if have_fcn('octave', 'vnum') >= 4
        file_in_path_warn_id = 'Octave:data-file-in-path';
    else
        file_in_path_warn_id = 'Octave:fopen-file-in-path';
    end
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

if ~have_fcn('regexp_split')
    t_skip(num_tests, 'PSSE2MPC requires newer MATLAB/Octave with regexp split support');
else
    t = '[records, sections] = psse_read() : length(records)';
    [records, sections] = psse_read(raw, verbose);
    t_is(length(records), 11, 12, t);
    t = '[records, sections] = psse_read() : length(sections)';
    t_is(length(sections), 3, 12, t);

    expected = { ...
        {1, 'Line 1   ', 1.1, -0.1, 0.011, 1, 1.1, 'A', '', 'A'}, ...
        {2, 'Line, "2"', 2.2, -0.2, 0.022, 2, 2.2, 'B', '', 'B'}, ...
        {3, 'Line, ''3''', 3.3, -0.3, 0.033, 3, 3.3, 'C', '', 'C'}, ...
        {4, sprintf('Line\t4'), 4.4, -0.4, 0.044, 4, 4.4, 'D', '', 'D'}, ...
    };
    ec = { ...
        ', "comment 1"', ...
        'comment, ''2''', ...
        sprintf('''comment\t3'''), ...
        '//comment,4', ...
    };

    for i = 1:sections(2).last - sections(2).first + 1
        t = sprintf('psse_parse_line(str%d, template) : ', i);
        [d, c] = psse_parse_line(records{i+sections(2).first-1}, 'dsffgDFcsc');
        t_is(length(d), length(expected{i}), 12, [t 'length']);
        for k = 1:length(d)
            if isnumeric(expected{i}{k})
                t_is(d{k}, expected{i}{k}, 12, sprintf('%s col %d', t, k));
            elseif isempty(expected{i}{k})
                t_ok(isempty(d{k}), sprintf('%s col %d', t, k));
            else
                t_ok(strcmp(d{k}, expected{i}{k}), sprintf('%s col %d', t, k));
            end
        end
        t_ok(strcmp(c, ec{i}), sprintf('%s comment', t));
    end

    t = 'psse_parse_line : missing optional columns : ';
    [d, c] = psse_parse_line(records{1}, 'dsffgDFcscdfgcs');
    t_is(length(d), 15, 12, [t 'length']);
    t_ok(all(cellfun(@isempty, d(11:15))), [t 'all empty']);

    t = 'psse_parse_section : ';
    [d, w] = psse_parse_section({}, records, sections, 2, 0, 'test1', 'dsFfgDF.sc');
    t_ok(isstruct(d) && isfield(d, 'num') && isfield(d, 'txt'), [t 'struct']);
    t_is(size(d.num), [4 11], 12, [t 'size(num)']);
    t_is(size(d.txt), [4 11], 12, [t 'size(txt)']);
    for i = 1:size(d.num, 1)
        for k = 1:size(d.num, 2)-1
            if isnumeric(expected{i}{k})
                t_is(d.num(i,k), expected{i}{k}, 12, sprintf('%s num(%d,%d)', t, i, k));
                t_ok(isempty(d.txt{i,k}), sprintf('%s txt{%d,%d}', t, i, k));
            elseif isempty(expected{i}{k})
                t_ok(isnan(d.num(i,k)), sprintf('%s num(%d,%d)', t, i, k));
                t_ok(isempty(d.txt{i,k}), sprintf('%s txt{%d,%d}', t, i, k));
            else
                t_ok(isnan(d.num(i,k)), sprintf('%s num(%d,%d)', t, i, k));
                t_ok(strcmp(d.txt{i,k}, expected{i}{k}), sprintf('%s txt{%d,%d}', t, i, k));
            end
        end
    end

    t = 'psse2mpc(rawfile, casefile)';
    txt = 'MATPOWER 5.0 using PSSE2MPC on 11-Aug-2014';
    for k = 2:3
        fname = sprintf(case_n, k);
        rawname = sprintf('%s.raw', fname);
        casename = sprintf('%s.m', fname);
        tmpfname = sprintf('%s_%d', fname, fix(1e9*rand));
        tmpcasename = sprintf('%s.m', tmpfname);
        mpc = psse2mpc(rawname, tmpfname, 0);
        str = fileread(casename);
        str2 = fileread(tmpcasename);
        str = strrep(str, char([13 10]), char(10));     %% Win to Unix EOL chars
        str2 = strrep(str2, char([13 10]), char(10));   %% Win to Unix EOL chars
        str2 = strrep(str2, 'e-005', 'e-05');           %% needed on Windoze, who knows why?
        str2 = strrep(str2, tmpfname, fname);
        str2 = strrep(str2, upper(tmpfname), upper(fname));
        str2 = regexprep(str2, 'MATPOWER (.*) using PSSE2MPC on \d\d-...-\d\d\d\d', txt);
        delete(tmpcasename);
        t_ok(strcmp(str, str2), sprintf('%s : %s', t, fname));
    end

    t = 'save2psse -> psse2mpc : ';
    mpc0 = loadcase(casefile);
    tmpfname = sprintf('t_save2psse_%d', fix(1e9*rand));
    tmpfname = save2psse(tmpfname, mpc0);
    mpc = psse2mpc(tmpfname, 0);
    delete(tmpfname);
    t_is(mpc.bus, mpc0.bus, 12, [t 'bus']);
    t_is(mpc.branch, mpc0.branch, 12, [t 'branch']);
    t_is(mpc.gen, mpc0.gen, 12, [t 'gen']);
    t_is(mpc.dcline, mpc0.dcline, 4, [t 'dcline']);
    t_ok(isequal(mpc.bus_name, mpc0.bus_name), [t 'bus_name']);
end

if have_fcn('octave')
    warning(s1.state, file_in_path_warn_id);
end

t_end;
