function test_matpower(verbose)
%TEST_MATPOWER  Run all MATPOWER tests.
%   test_matpower runs all of the MATPOWER tests.
%   test_matpower(verbose) prints the details of the individual tests
%   if verbose is true.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_loadcase';
tests{end+1} = 't_ext2int2ext';
tests{end+1} = 't_jacobian';
tests{end+1} = 't_hessian';
tests{end+1} = 't_hasPQcap';
tests{end+1} = 't_pf';
if have_fcn('fmincon')
    tests{end+1} = 't_opf_fmincon';
end
if have_fcn('minopf')
    tests{end+1} = 't_opf_minopf';
end
tests{end+1} = 't_opf_pdipm';
tests{end+1} = 't_opf_scpdipm';
if have_fcn('pdipmopf')
    tests{end+1} = 't_opf_tspopf_pdipm';
end
if have_fcn('scpdipmopf')
    tests{end+1} = 't_opf_tspopf_scpdipm';
end
if have_fcn('tralmopf')
    tests{end+1} = 't_opf_tspopf_tralm';
if have_fcn('constr')
    tests{end+1} = 't_opf_constr';
end
if have_fcn('bpmpd')
    tests{end+1} = 't_opf_lp_den';
    tests{end+1} = 't_opf_lp_spr';
    tests{end+1} = 't_opf_lp_spf';
    tests{end+1} = 't_opf_dc_bpmpd';
end
if have_fcn('quadprog')
    tests{end+1} = 't_opf_dc_ot';
end
tests{end+1} = 't_opf_dc_pdipm';
tests{end+1} = 't_opf_dc_scpdipm';
tests{end+1} = 't_runopf_w_res';
tests{end+1} = 't_makePTDF';
tests{end+1} = 't_makeLODF';
tests{end+1} = 't_total_load';
tests{end+1} = 't_scale_load';

%% smartmarket tests
if exist('runmarket') == 2
    tests{end+1} = 't_off2case';
    if have_fcn('minopf')
        tests{end+1} = 't_auction_minopf';
    end
    tests{end+1} = 't_auction_pdipm';
    if have_fcn('pdipmopf')
        tests{end+1} = 't_auction_tspopf_pdipm';
    end
    tests{end+1} = 't_runmarket';
end

%% sopf tests
if exist('apply_contingency') == 2
    tests{end+1} = 't_apply_contingency';
end

t_run_tests( tests, verbose );

return;
