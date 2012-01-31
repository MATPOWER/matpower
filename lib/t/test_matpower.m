function test_matpower(verbose)
%TEST_MATPOWER  Run all MATPOWER tests.
%   TEST_MATPOWER runs all of the MATPOWER tests.
%   TEST_MATPOWER(VERBOSE) prints the details of the individual tests
%   if VERBOSE is true.
%
%   See also T_RUN_TESTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    verbose = 0;
end

tests = {};

%% MATPOWER base test
tests{end+1} = 't_loadcase';
tests{end+1} = 't_ext2int2ext';
tests{end+1} = 't_jacobian';
tests{end+1} = 't_hessian';
tests{end+1} = 't_totcost';
tests{end+1} = 't_modcost';
tests{end+1} = 't_hasPQcap';
if have_fcn('anon_fcns')
    tests{end+1} = 't_mips';
else
    tests{end+1} = 't_mips6';
end
tests{end+1} = 't_qps_matpower';
tests{end+1} = 't_pf';
tests{end+1} = 't_islands';
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
if have_fcn('constr')
    tests{end+1} = 't_opf_constr';
end
if have_fcn('bpmpd')
    tests{end+1} = 't_opf_lp_den';
    tests{end+1} = 't_opf_lp_spr';
    tests{end+1} = 't_opf_lp_spf';
    tests{end+1} = 't_opf_dc_bpmpd';
end
if have_fcn('cplex')
    tests{end+1} = 't_opf_dc_cplex';
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
tests{end+1} = 't_opf_userfcns';
tests{end+1} = 't_runopf_w_res';
tests{end+1} = 't_dcline';
tests{end+1} = 't_makePTDF';
tests{end+1} = 't_makeLODF';
tests{end+1} = 't_total_load';
tests{end+1} = 't_scale_load';

%% smartmarket tests
if exist('runmarket', 'file') == 2
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

%% sopf tests
if exist('apply_contingency', 'file') == 2
    tests{end+1} = 't_apply_contingency';
    tests{end+1} = 't_c3sopf_sopf2';
end

t_run_tests( tests, verbose );
