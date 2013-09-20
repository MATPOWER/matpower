function [cb_state, results] = ...
    cpf_default_callback(k, V_c, lam_c, V_p, lam_p, cb_data, cb_state, results)
%CPF_DEFAULT_CALLBACK   Default callback function for CPF
%   [CB_STATE, RESULTS] = ...
%       CPF_DEFAULT_CALLBACK(K, V_C, LAM_C, V_P, LAM_P, ...
%                            CB_DATA, CB_STATE, RESULTS)
%
%   Default callback function used by RUNCPF. Takes input from current
%   iteration, returns a user defined state struct and on the final call
%   a results struct.
%
%   Inputs:
%       K - continuation step iteration count
%       V_C - vector of complex bus voltages after K-th corrector step
%       LAM_C - value of LAMBDA after K-th corrector step
%       V_P - vector of complex bus voltages after K-th predictor step
%       LAM_P - value of LAMBDA after K-th predictor step
%       CB_DATA - struct containing potentially useful static data,
%           with the following fields (all based on internal indexing):
%           mpc_base - MATPOWER case struct of base state
%           mpc_target - MATPOWER case struct of target state
%           Ybus - bus admittance matrix
%           Yf - branch admittance matrix, "from" end of branches
%           Yt - branch admittance matrix, "to" end of branches
%           pv - list of indices of PV buses
%           pq - list of indices of PQ buses
%           ref - list of indices of REF buses
%           mpopt - MATPOWER options vector
%       CB_STATE - user-defined struct containing any
%           information the callback function would like to
%           pass from one invokation to the next
%       RESULTS - initial value of output struct to be assigned to
%           CPF field of results struct returned by RUNCPF
%
%   Outputs:
%       CB_STATE - updated version of CB_STATE input arg
%       RESULTS - updated version of RESULTS input arg
%
%   This function is called in three different contexts, distinguished
%   as follows:
%   (1) initial - called without RESULTS output arg, with K = 0,
%           after base power flow, before 1st CPF step.
%   (2) iterations - called without RESULTS output arg, with K > 0
%           at each iteration, after predictor-corrector step
%   (3) final - called with RESULTS output arg, after exiting
%           predictor-corrector loop, inputs identical to
%           last iteration call
%
%   See also RUNCPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2013 by Power System Engineering Research Center (PSERC)
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

%%-----  INITIAL call  -----
if k == 0
    cb_state = struct(  'V_p', V_p, ...
                        'lam_p', lam_p, ...
                        'V_c', V_c, ...
                        'lam_c', lam_c, ...
                        'iterations', 0);
%%-----  FINAL call  -----
elseif nargout == 2
    results = cb_state; %% initialize results with final state
    results.max_lam = max(cb_state.lam_c);
    results.iterations = k;
%%-----  ITERATION call  -----
else
    cb_state.V_p   = [cb_state.V_p V_p];
    cb_state.lam_p = [cb_state.lam_p lam_p];
    cb_state.V_c   = [cb_state.V_c V_c];
    cb_state.lam_c = [cb_state.lam_c lam_c];
    cb_state.iterations    = k;
end
