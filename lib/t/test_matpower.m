function test_matpower(verbose)
%TEST_MATPOWER  Run all MATPOWER tests.
%   TEST_MATPOWER runs all of the MATPOWER tests.
%   TEST_MATPOWER(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_test_fcns';
tests{end+1} = 't_nested_struct_copy';
tests{end+1} = 't_feval_w_path';
tests{end+1} = 't_mpoption';
tests{end+1} = 't_loadcase';
tests{end+1} = 't_ext2int2ext';
tests{end+1} = 't_jacobian';
tests{end+1} = 't_hessian';
tests{end+1} = 't_margcost';
tests{end+1} = 't_totcost';
tests{end+1} = 't_modcost';
tests{end+1} = 't_hasPQcap';
tests{end+1} = 't_mplinsolve';
tests{end+1} = 't_mips';
tests{end+1} = 't_qps_matpower';
tests{end+1} = 't_miqps_matpower';
tests{end+1} = 't_pf';
tests{end+1} = 't_cpf';
tests{end+1} = 't_islands';
tests{end+1} = 't_opf_model';
if have_fcn('fmincon')
    tests{end+1} = 't_opf_fmincon';
end
if have_fcn('minopf')
    tests{end+1} = 't_opf_minopf';
end
tests{end+1} = 't_opf_mips';
tests{end+1} = 't_opf_mips_sc';
if have_fcn('pdipmopf')
    tests{end+1} = 't_opf_tspopf_pdipm';
end
if have_fcn('scpdipmopf')
    tests{end+1} = 't_opf_tspopf_scpdipm';
end
if have_fcn('tralmopf')
    tests{end+1} = 't_opf_tspopf_tralm';
end
if have_fcn('ipopt')
    tests{end+1} = 't_opf_ipopt';
end
if have_fcn('knitro')
    tests{end+1} = 't_opf_knitro';
end
if have_fcn('bpmpd')
    tests{end+1} = 't_opf_dc_bpmpd';
end
if have_fcn('clp')
    tests{end+1} = 't_opf_dc_clp';
end
if have_fcn('cplex')
    tests{end+1} = 't_opf_dc_cplex';
end
if have_fcn('glpk')
    tests{end+1} = 't_opf_dc_glpk';
end
if have_fcn('gurobi')
    tests{end+1} = 't_opf_dc_gurobi';
end
if have_fcn('ipopt')
    tests{end+1} = 't_opf_dc_ipopt';
end
tests{end+1} = 't_opf_dc_mips';
tests{end+1} = 't_opf_dc_mips_sc';
if have_fcn('mosek')
    tests{end+1} = 't_opf_dc_mosek';
end
if have_fcn('quadprog')
    tests{end+1} = 't_opf_dc_ot';
end
if have_fcn('sdp_pf')
    tests{end+1} = 't_opf_sdpopf';
    tests{end+1} = 't_testglobalopt';
    tests{end+1} = 't_insolvablepf';
    tests{end+1} = 't_insolvablepf_limitQ';
    tests{end+1} = 't_insolvablepfsos';
    tests{end+1} = 't_insolvablepfsos_limitQ';
end
tests{end+1} = 't_opf_userfcns';
tests{end+1} = 't_opf_softlims';
tests{end+1} = 't_runopf_w_res';
tests{end+1} = 't_dcline';
tests{end+1} = 't_get_losses';
tests{end+1} = 't_makePTDF';
tests{end+1} = 't_makeLODF';
tests{end+1} = 't_printpf';
tests{end+1} = 't_vdep_load';
tests{end+1} = 't_total_load';
tests{end+1} = 't_scale_load';
tests{end+1} = 't_apply_changes';
tests{end+1} = 't_psse';

%% smartmarket tests
if have_fcn('smartmarket')
    tests{end+1} = 't_off2case';
    if have_fcn('minopf')
        tests{end+1} = 't_auction_minopf';
    end
    tests{end+1} = 't_auction_mips';
    if have_fcn('pdipmopf')
        tests{end+1} = 't_auction_tspopf_pdipm';
    end
    tests{end+1} = 't_runmarket';
end

t_run_tests( tests, verbose );
