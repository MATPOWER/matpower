function t_run_tests(test_names, verbose)
%T_RUN_TESTS  Run a series of tests.
%   t_run_tests(test_names, verbose) runs a set of tests whose names
%   are given in the cell array test_names. If the optional parameter
%   verbose is true, it prints the details of the individual tests.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

if nargin < 2
	verbose = 0;
end

global t_num_of_tests;
global t_counter;
global t_ok_cnt;
global t_not_ok_cnt;

%% figure out padding for printing
if ~verbose
	maxlen = 0;
	len = zeros(length(test_names), 1);
	for k = 1:length(test_names)
		len(k) = length(test_names{k});
	end
	maxlen = max(len);
end

%% initialize statistics
num_of_tests = 0;
counter = 0;
ok_cnt = 0;
not_ok_cnt = 0;

t0 = clock;
for k = 1:length(test_names)
	if verbose
		fprintf('\n----------  %s  ----------\n', test_names{k});
	else
		pad = maxlen + 4 - length(test_names{k});
		fprintf('%s', test_names{k});
		for m = 1:pad, fprintf('.'); end
	end
	feval( test_names{k}, ~verbose );
	
	num_of_tests	= num_of_tests	+ t_num_of_tests;
	counter			= counter		+ t_counter;
	ok_cnt			= ok_cnt		+ t_ok_cnt;
	not_ok_cnt		= not_ok_cnt	+ t_not_ok_cnt;
end

if verbose
	fprintf('\n\n----------  Summary  ----------\n');
end
if counter == num_of_tests & ok_cnt == counter & not_ok_cnt == 0
	fprintf('All tests successful (%d of %d).\n', ok_cnt, num_of_tests);
else
	fprintf('Ran %d of %d tests: %d passed, %d failed.\n', ...
		counter, num_of_tests, ok_cnt, not_ok_cnt);
end
fprintf('Elapsed time %.2f seconds.\n', etime(clock, t0));


return;
