function test_matpower(verbose)
%TEST_ALL  Run all MATPOWER tests.
%   test_matpower runs all of the MATPOWER tests.
%   test_matpower(verbose) prints the details of the individual tests
%   if verbose is true.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

if nargin < 1
	verbose = 0;
end

test_list = {	't_loadcase', ...
				't_jacobian', ...
				't_pf', ...
				't_opf'	};

%% add smartmarket tests if available
nt = size(test_list, 2);
if exist('runmkt')
	test_list{1,nt+1} = 't_auction';
end

t_run_tests( test_list, verbose );

return;
