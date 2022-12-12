function ok = t_file_match(got_fname, exp_fname, msg, reps, del_got_fname)
%T_FILE_MATCH  Tests whether two text files match.
%   T_FILE_MATCH(GOT_FNAME, EXP_FNAME, MSG)
%   T_FILE_MATCH(GOT_FNAME, EXP_FNAME, MSG, REPS)
%   T_FILE_MATCH(GOT_FNAME, EXP_FNAME, MSG, REPS, DEL_GOT_FNAME)
%   Uses T_STR_MATCH on the contents of two text files whos names/paths
%   are given in GOT_FNAME and EXP_FNAME. If both files exist and the
%   contents match, the test passes.
%
%   It ignores any differences in line ending characters and, like
%   T_STR_MATCH, can apply replacements to the contents of GOT_FNAME,
%   and optionally EXP_FNAME, as specified by REPS before comparing.
%
%   The REPS argument is a cell array of replacement specs, applied
%   sequentially, where each replacement spec is a cell array of the
%   following form:
%       {ORIGINAL, REPLACEMENT}
%       {ORIGINAL, REPLACEMENT, RE}
%       {ORIGINAL, REPLACEMENT, RE, BOTH}
%   The ORIGINAL and REPLACEMENT arguments are passed directly as the
%   2nd and 3rd arguments to REGEXPREP (or to STRREP if RE is present
%   and false). The replacement applies to GOT_FNAME only, unless BOTH is
%   present and true, in which case it also applies to EXP_FNAME.
%
%   If DEL_GOT_FNAME is present and true it will delete the file named
%   in GOT_FNAME if the test passes.
%
%   Optionally returns a true or false value indicating whether or
%   not the test succeeded.
%
%   Example:
%       quiet = 0;
%       t_begin(5, quiet);
%
%       % replace Windows EOL chars with Unix EOL chars
%       reps = {{char([13 10]), char(10), 0, 1}};
%       got_fname = 'test_generated_output.txt';
%       exp_fname = 'expected_output.txt';
%       t_file_match(got_fname, exp_fname, 'mytest', reps, 1);
%
%       t_end;
%
%   See also T_OK, T_IS, T_STR_MATCH, T_SKIP, T_BEGIN, T_END, T_RUN_TESTS.

%   MP-Test
%   Copyright (c) 2004-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

%% default arguments
if nargin < 5
    del_got_fname = 0;
    if nargin < 4
        reps = {};
    end
end

%% check existence
got_exists = exist(got_fname) == 2;

%% compare contents
if got_exists
    %% read both files
    got = fileread(got_fname);
    expected = fileread(exp_fname);

    %% ignore line endings by replacing Win EOL with Unix EOL chars
    reps = {{char([13 10]), char(10), 0, 1}, reps{:}};

    %% check if contents match
    TorF = t_str_match(got, expected, msg, reps);
    
    %% delete GOT_FNAME if requested and test passed
    if TorF && del_got_fname
        delete(got_fname);
    end
else
    TorF = t_ok(got_exists, [msg ' - missing file']);
end

if nargout
    ok = TorF;
end
