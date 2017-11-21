function qpopt = mpopt2qpopt(mpopt, model, alg)
%MPOPT2QPOPT   Create/modify MI/QPS_MATPOWER options struct from MPOPT.
%
%   QPOPT = MPOPT2QPOPT(MPOPT, MODEL)
%   QPOPT = MPOPT2QPOPT(MPOPT, MODEL, ALG)
%
%   Uses a MATPOWER options struct, MPOPT, to create or modify an
%   MIQPS_MATPOWER or QPS_MATPOWER options struct.
%
%   Inputs (default values in parentheses):
%       MPOPT : MATPOWER options struct
%       MODEL ('MIQP') : (optional) one of the following model types, required
%               for selection of solver in case ALG is 'DEFAULT' (solver
%               precedence for each model type list in parentheses):
%           'LP'   - linear program with all continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT (if MATLAB), GLPK, BPMPD, MIPS)
%           'QP'   - quadratic program with all continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT (if large-scale alg available),
%                    BPMPD, MIPS)
%           'MILP' - LP with mixed integer/continuous variables
%                   (GUROBI, CPLEX, MOSEK, OT, GLPK)
%           'MIQP' - (default) QP with mixed integer/continuous variables
%                   (GUROBI, CPLEX, MOSEK)
%       ALG ('opf.dc') : (optional) 'opf.dc', 'most', or any valid value of
%               OPT.alg for QPS_MATPOWER or MIQPS_MATPOWER. The first two
%               options indicate that it should be taken from
%               MPOPT.opf.dc.solver or MPOPT.most.solver, respectively.
%
%   Output:
%       QPOPT : an options struct for use by QPS_MATPOWER or MIQPS_MATPOWER
%               and friends

%   MATPOWER
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% set default args
if nargin < 3
    alg = '';
    if nargin < 2
        model = '';
    end
end
if isempty(model)
    model = 'MIQP';
else
    model = upper(model);
end
skip_prices = 0;
price_stage_warn_tol = [];

%% get ALG from mpopt, if necessary
switch alg
    case {'opf.dc', ''}
        alg = upper(mpopt.opf.dc.solver);
    case 'most'
        alg = upper(mpopt.most.solver);
        skip_prices             = mpopt.most.skip_prices;
        price_stage_warn_tol    = mpopt.most.price_stage_warn_tol;
    otherwise
        alg = upper(alg);
end

%% default solver
switch alg
    case {'DEFAULT', 0}
        if have_fcn('gurobi')
            alg = 'GUROBI';     %% use Gurobi by default, if available
        elseif have_fcn('cplex')
            alg = 'CPLEX';      %% if not, then CPLEX, if available
        elseif have_fcn('mosek')
            alg = 'MOSEK';      %% if not, then MOSEK, if available
        elseif have_fcn('linprog_ds') && strcmp(model, 'LP') && have_fcn('matlab') || ...
                have_fcn('quadprog_ls') && strcmp(model, 'QP') || ...
                have_fcn('intlinprog') && strcmp(model, 'MILP')
            alg = 'OT';         %% if not, then newer Optimization Tbx, if
                                %% available and applicable
        elseif have_fcn('glpk') && model(end-1) == 'L'  %% LP or MILP
            alg = 'GLPK';       %% if not, then GLPK, if available & applicable
        elseif have_fcn('linprog') && strcmp(model, 'LP') && have_fcn('matlab')
            alg = 'OT';         %% if not, then older Optimization Tbx, if
                                %% available and applicable
        elseif model(1) ~= 'M'  %% LP or QP
            if have_fcn('bpmpd')
                alg = 'BPMPD';  %% if not, then BPMPD_MEX, if available
                                %% and applicable
            else
                alg = 'MIPS';   %% otherwise MIPS, if applicable
            end
        else
            error('mpopt2qpopt: Sorry, no solver available for %s models', model);
        end
end

%% create MI/QPS_MATPOWER options struct
qpopt = struct('alg', alg, 'verbose', mpopt.verbose);
switch alg
    case {'MIPS', 200, 250}
        %% set up options
        qpopt.mips_opt = mpopt.mips;
        if qpopt.mips_opt.feastol == 0      %% = MPOPT.opf.violation by default
            qpopt.mips_opt.feastol = mpopt.opf.violation;
        end
    case {'IPOPT', 400}
        qpopt.ipopt_opt = ipopt_options([], mpopt);
    case {'BPMPD', 100}
        bp_opt = bpopt;
        bp_opt(20) = 1e-9;  %% TOPT1
        qpopt.bp_opt = bp_opt;
    case 'CLP'
        qpopt.clp_opt = clp_options([], mpopt);
    case {'CPLEX', 500}
        qpopt.cplex_opt = cplex_options([], mpopt);
    case 'GLPK'
        qpopt.glpk_opt = glpk_options([], mpopt);
    case {'GUROBI', 700}
        qpopt.grb_opt = gurobi_options([], mpopt);
    case {'MOSEK', 600}
        qpopt.mosek_opt = mosek_options([], mpopt);
    case {'OT', 300}
        if isfield(mpopt, 'linprog') && ~isempty(mpopt.linprog)
            qpopt.linprog_opt = mpopt.linprog;
        end
        if isfield(mpopt, 'quadprog') && ~isempty(mpopt.quadprog)
            qpopt.quadprog_opt = mpopt.quadprog;
        end
        if isfield(mpopt, 'intlinprog') && ~isempty(mpopt.intlinprog)
            qpopt.intlinprog_opt = mpopt.intlinprog;
        end
end
if model(1) == 'M'
    qpopt.skip_prices           = skip_prices;
    qpopt.price_stage_warn_tol  = price_stage_warn_tol;
end
