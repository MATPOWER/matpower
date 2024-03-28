function nleqopt = mpopt2nleqopt(mpopt, model, alg)
% mpopt2nleqopt - Create/modify nleqs_master options struct from ``mpopt``.
% ::
%
%   NLEQOPT = MPOPT2NLEQOPT(MPOPT, MODEL)
%   NLEQOPT = MPOPT2NLEQOPT(MPOPT, MODEL, ALG)
%
%   Uses a MATPOWER options struct, MPOPT, to create or modify an
%   NLEQS_MASTER options struct.
%
%   Inputs (default values in parentheses):
%       MPOPT : MATPOWER options struct
%       MODEL ('NLEQ') : (optional) one of the following model types, required
%               for selection of solver in case ALG is 'DEFAULT' (solver
%               precedence for each model type list in parentheses):
%           'NLEQ' - nonlinear equation with all continuous variables
%                   (NEWTON, FSOLVE)
%       ALG ('pf.alg') : (optional) 'pf.alg' or any valid value of
%               OPT.alg for NLEQS_MASTER. The first option indicates that
%               it should be taken from MPOPT.pf.alg.
%
%   Output:
%       NLEQOPT : an options struct for use by NLEQS_MASTER and friends
%
% See also nleqs_master, mpoption.

%   MP-Opt-Model
%   Copyright (c) 2015-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% set default args
if nargin < 3
    alg = '';
    if nargin < 2
        model = '';
    end
end
if isempty(model)
    model = 'NLEQ';
else
    model = upper(model);
end

%% get ALG, MAX_IT from mpopt, if necessary
switch alg
    case {'pf.alg', ''}
        alg = upper(mpopt.pf.alg);
    otherwise
        alg = upper(alg);
end

%% create NLEQS_MASTER options struct
nleqopt = struct(   'verbose',  mpopt.verbose, ...
                    'alg',      alg, ...
                    'tol',      mpopt.pf.tol  );
switch alg
    case {'DEFAULT', 'NEWTON'}
        %% set up options
        nleqopt.newton_opt.lin_solver = mpopt.pf.nr.lin_solver;
        nleqopt.max_it = mpopt.pf.nr.max_it;
    case 'FD'
        nleqopt.max_it = mpopt.pf.fd.max_it;
    case 'GS'
        nleqopt.max_it = mpopt.pf.gs.max_it;
    case 'ZG'
        nleqopt.max_it = mpopt.pf.zg.max_it;
        nleqopt.alg = 'CORE';
    case 'FSOLVE'
        %% basic optimset options needed for fmincon
%         nleqopt.fsolve_opt.Algorithm = '';
%         nleqopt.fsolve_opt.Algorithm = 'trust-region-dogleg';
%         nleqopt.fsolve_opt.Algorithm = 'trust-region';              %% for optimoptions
%         nleqopt.fsolve_opt.Algorithm = 'trust-region-reflective';   %% for optimset
%         nleqopt.fsolve_opt.Algorithm = 'levenberg-marquardt';
        nleqopt.max_it = 0;      %% use fsolve default
end
