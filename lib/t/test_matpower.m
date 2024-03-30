function success = test_matpower(verbose, exit_on_fail)
% test_matpower - Run all |MATPOWER| tests.
% ::
%
%   test_matpower
%   test_matpower(verbose)
%   test_matpower(verbose, exit_on_fail)
%   success = test_matpower(...)
%
% Runs all of the |MATPOWER| tests. If ``verbose`` is true *(false by
% default)*, it prints the details of the individual tests. If
% ``exit_on_fail`` is true *(false by default)*, it will exit MATLAB or
% Octave with a status of 1 unless t_run_tests returns ``all_ok`` true.
%
% See also t_run_tests.

%   MATPOWER
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 2
    exit_on_fail = 0;
    if nargin < 1
        verbose = 0;
    end
end

tests = {};

%% MP-Test
tests{end+1} = 't_test_fcns';
tests{end+1} = 't_have_feature';

%% MIPS
tests{end+1} = 't_mplinsolve';
tests{end+1} = 't_mips';
tests{end+1} = 't_mips_pardiso';
tests{end+1} = 't_qps_mips';

%% MP-Opt-Model
tests{end+1} = 't_have_fcn';
tests{end+1} = 't_nested_struct_copy';
tests{end+1} = 't_nleqs_master';
tests{end+1} = 't_pnes_master';
tests{end+1} = 't_qps_master';
tests{end+1} = 't_miqps_master';
tests{end+1} = 't_nlps_master';
tests{end+1} = 't_opt_model';
tests{end+1} = 't_om_solve_leqs';
tests{end+1} = 't_om_solve_nleqs';
tests{end+1} = 't_om_solve_pne';
tests{end+1} = 't_om_solve_qps';
tests{end+1} = 't_om_solve_miqps';
tests{end+1} = 't_om_solve_nlps';

%% deprecated wrappers
tests{end+1} = 't_qps_matpower';
tests{end+1} = 't_miqps_matpower';

%% MATPOWER
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

%% MP-Core
if have_feature('mp_core')
    tests{end+1} = 't_mp_mapped_array';
    tests{end+1} = 't_mp_table';
    tests{end+1} = 't_mp_data_model';
    tests{end+1} = 't_dmc_element';
    tests{end+1} = 't_mp_dm_converter_mpc2';
    tests{end+1} = 't_nm_element';
    tests{end+1} = 't_port_inj_current_acc';
    tests{end+1} = 't_port_inj_current_acp';
    tests{end+1} = 't_port_inj_power_acc';
    tests{end+1} = 't_port_inj_power_acp';
    tests{end+1} = 't_node_test';
    tests{end+1} = 't_run_mp';
    tests{end+1} = 't_run_opf_default';
    tests{end+1} = 't_run_mp_3p';
    tests{end+1} = 't_pretty_print';
    tests{end+1} = 't_mpxt_reserves';
    tests{end+1} = 't_mpxt_legacy_dcline';
end

%% more MATPOWER
tests{end+1} = 't_pf_ac';
tests{end+1} = 't_pf_dc';
tests{end+1} = 't_pf_radial';
tests{end+1} = 't_cpf';
tests{end+1} = 't_islands';
tests{end+1} = 't_opf_model';
tests{end+1} = 't_opf_default';
if have_feature('fmincon')
    tests{end+1} = 't_opf_fmincon';
end
if have_feature('minopf')
    tests{end+1} = 't_opf_minopf';
end
tests{end+1} = 't_opf_mips';
if have_feature('pdipmopf')
    tests{end+1} = 't_opf_tspopf_pdipm';
end
if have_feature('scpdipmopf')
    tests{end+1} = 't_opf_tspopf_scpdipm';
end
if have_feature('tralmopf')
    tests{end+1} = 't_opf_tspopf_tralm';
end
if have_feature('ipopt')
    tests{end+1} = 't_opf_ipopt';
end
if have_feature('knitro')
    tests{end+1} = 't_opf_knitro';
end
tests{end+1} = 't_opf_dc_default';
if have_feature('bpmpd')
    tests{end+1} = 't_opf_dc_bpmpd';
end
if have_feature('clp')
    tests{end+1} = 't_opf_dc_clp';
end
if have_feature('cplex')
    tests{end+1} = 't_opf_dc_cplex';
end
if have_feature('glpk')
    tests{end+1} = 't_opf_dc_glpk';
end
if have_feature('gurobi')
    tests{end+1} = 't_opf_dc_gurobi';
end
if have_feature('ipopt')
    tests{end+1} = 't_opf_dc_ipopt';
end
tests{end+1} = 't_opf_dc_mips';
tests{end+1} = 't_opf_dc_mips_sc';
if have_feature('mosek')
    tests{end+1} = 't_opf_dc_mosek';
end
if have_feature('osqp')
    tests{end+1} = 't_opf_dc_osqp';
end
if have_feature('quadprog')
    tests{end+1} = 't_opf_dc_ot';
end

%% SDP_PF
if have_feature('sdp_pf')
    if have_feature('mosek') || have_feature('sdpt3') || have_feature('sedumi')
        tests{end+1} = 't_opf_sdpopf';
        tests{end+1} = 't_insolvablepf';
        tests{end+1} = 't_insolvablepf_limitQ';
        tests{end+1} = 't_insolvablepfsos';
        tests{end+1} = 't_insolvablepfsos_limitQ';
    end
    tests{end+1} = 't_testglobalopt';
end
tests{end+1} = 't_opf_userfcns';
tests{end+1} = 't_opf_softlims';
tests{end+1} = 't_runopf_w_res';
tests{end+1} = 't_dcline';
tests{end+1} = 't_get_losses';
tests{end+1} = 't_load2disp';
tests{end+1} = 't_makePTDF';
tests{end+1} = 't_makeLODF';
tests{end+1} = 't_printpf';
tests{end+1} = 't_vdep_load';
tests{end+1} = 't_total_load';
tests{end+1} = 't_scale_load';
tests{end+1} = 't_apply_changes';
tests{end+1} = 't_psse';

%% smartmarket tests
if have_feature('smartmarket')
    tests{end+1} = 't_off2case';
    if have_feature('minopf')
        tests{end+1} = 't_auction_minopf';
    end
    tests{end+1} = 't_auction_mips';
    if have_feature('pdipmopf')
        tests{end+1} = 't_auction_tspopf_pdipm';
    end
    tests{end+1} = 't_runmarket';
end

%% run the tests
all_ok = t_run_tests( tests, verbose );

%% handle success/failure
if exit_on_fail && ~all_ok
    exit(1);
end
if nargout
    success = all_ok;
end
