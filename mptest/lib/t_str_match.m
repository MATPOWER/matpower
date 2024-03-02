function ok = t_str_match(got, expected, msg, reps)
% t_str_match - Tests whether two strings/char arrays match.
% ::
%
%   t_str_match(got, expected, msg)
%   t_str_match(got, expected, msg, reps)
%   ok = t_str_match(...)
%
% This is equivalent to ``t_ok(strcmp(got, expected), msg)``, with the
% option to apply replacements to ``got``, and optionally to ``expected``,
% as specified by ``reps`` before comparing.
%
% Inputs:
%   got (char array) : actual result
%   expected (char array) : expected result
%   msg (char array) : message to display for this test
%   reps (char array) : *(optional)* cell array of replacement specs
%
% Output:
%   ok (boolean) : *(optional)* true if test passed, false if failed
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
% is present and false). The replacement applies to ``got`` only,
% unless ``both`` is present and true, in which case it also applies to
% ``expected``.
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
%   got = fileread('test_generated_output.txt');
%   expected = fileread('expected_output.txt');
%   t_str_match(got, expected, 'mytest', reps);
%
%   t_end;
%
% See also t_ok, t_is, t_file_match, t_skip, t_begin, t_end, t_run_tests.

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

%% default arguments
if nargin < 4
    reps = {};
end

%% do the replacments
for k = 1:length(reps)
    rep = reps{k};
    re = length(rep) < 3 || rep{3};     %% use regexprep()
    both = length(rep) > 3 && rep{4};   %% repl in expected as well as got
    if re       %% use regexprep()
        got = regexprep(got, rep{1}, rep{2});
        if both
            expected = regexprep(expected, rep{1}, rep{2});
        end
    else        %% use strrep()
        got = strrep(got, rep{1}, rep{2});
        if both
            expected = strrep(expected, rep{1}, rep{2});
        end
    end
end

%% check for match
TorF = t_ok(strcmp(got, expected), msg);

if nargout
    ok = TorF;
end
