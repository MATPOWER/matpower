function ok = t_file_match(got_fname, exp_fname, msg, reps, del_got_fname)
% t_file_match - Tests whether two text files match.
% ::
%
%   t_file_match(got_fname, exp_fname, msg)
%   t_file_match(got_fname, exp_fname, msg, reps)
%   t_file_match(got_fname, exp_fname, msg, reps, del_got_fname)
%   ok = t_file_match(...)
%
% Uses t_str_match on the contents of two text files whos names/paths
% are given in ``got_fname`` and ``exp_fname``. If both files exist and the
% contents match, the test passes.
%
% Inputs:
%   got_fname (char array) : name of file containing actual result
%   exp_fname (char array) : name of file containing expected result
%   msg (char array) : message to display for this test
%   reps (char array) : *(optional)* cell array of replacement specs
%   del_got_fname (boolean) : *(optional, default = false)* if true, the
%       file named in ``got_fname`` will be deleted if the test passes
%
% Output:
%   ok (boolean) : *(optional)* true if test passed, false if failed
%
% It ignores any differences in line ending characters and, like
% t_str_match, can apply replacements to the contents of ``got_fname``,
% and optionally ``exp_fname``, as specified by ``reps`` before comparing.
%
% The ``reps`` argument is a cell array of replacement specs, applied
% sequentially, where each replacement spec is a cell array of the
% following form::
%
%   {original, replacement}
%   {original, replacement, re}
%   {original, replacement, re, both}
%
% The ``original`` and ``replacement`` arguments are passed directly as the
% 2nd and 3rd arguments to :func:`regexprep` (or to :func:`strrep` if ``re``
% is present and false). The replacement applies to ``got_fname`` only,
% unless ``both`` is present and true, in which case it also applies to
% ``exp_fname``.
%
% If ``del_got_fname`` is present and true it will delete the file named
% in ``got_fname`` if the test passes.
%
% Intended to be called between calls to t_begin and t_end.
%
% Example::
%
%   quiet = 0;
%   t_begin(5, quiet);
%
%   % replace Windows EOL chars with Unix EOL chars
%   reps = {{char([13 10]), char(10), 0, 1}};
%   got_fname = 'test_generated_output.txt';
%   exp_fname = 'expected_output.txt';
%   t_file_match(got_fname, exp_fname, 'mytest', reps, 1);
%
%   t_end;
%
% See also t_ok, t_is, t_str_match, t_skip, t_begin, t_end, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
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
