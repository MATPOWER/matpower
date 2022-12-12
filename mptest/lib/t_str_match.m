function ok = t_str_match(got, expected, msg, reps)
%T_STR_MATCH  Tests whether two strings/char arrays match.
%   T_STR_MATCH(GOT, EXPECTED, MSG)
%   T_STR_MATCH(GOT, EXPECTED, MSG, REPS)
%   This is equivalent to T_OK(STRCMP(GOT, EXPECTED), MSG), with the
%   option to apply replacements to GOT, and optionally to EXPECTED,
%   as specified by REPS before comparing.
%
%   The REPS argument is a cell array of replacement specs, applied
%   sequentially, where each replacement spec is a cell array of the
%   following form:
%       {ORIGINAL, REPLACEMENT}
%       {ORIGINAL, REPLACEMENT, RE}
%       {ORIGINAL, REPLACEMENT, RE, BOTH}
%   The ORIGINAL and REPLACEMENT arguments are passed directly as the
%   2nd and 3rd arguments to REGEXPREP (or to STRREP if RE is present
%   and false). The replacement applies to GOT only, unless BOTH is
%   present and true, in which case it also applies to EXPECTED.
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
%       got = fileread('test_generated_output.txt');
%       expected = fileread('expected_output.txt');
%       t_str_match(got, expected, 'mytest', reps);
%
%       t_end;
%
%   See also T_OK, T_IS, T_FILE_MATCH, T_SKIP, T_BEGIN, T_END, T_RUN_TESTS.

%   MP-Test
%   Copyright (c) 2004-2022, Power Systems Engineering Research Center (PSERC)
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
