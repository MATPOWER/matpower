function nlpopt = mpopt2nlpopt(mpopt, model, alg)
% mpopt2nlpopt - Create/modify nlps_master options struct from ``mpopt``.
% ::
%
%   NLPOPT = MPOPT2NLPOPT(MPOPT, MODEL)
%   NLPOPT = MPOPT2NLPOPT(MPOPT, MODEL, ALG)
%
%   Uses a MATPOWER options struct, MPOPT, to create or modify an
%   NLPS_MASTER options struct.
%
%   Inputs (default values in parentheses):
%       MPOPT : MATPOWER options struct
%       MODEL ('NLP') : (optional) one of the following model types, required
%               for selection of solver in case ALG is 'DEFAULT' (solver
%               precedence for each model type list in parentheses):
%           'NLP'   - nonlinear program with all continuous variables
%                   (MIPS, FMINCON, IPOPT, Artelys Knitro)
%           'MINLP' - NLP with mixed integer/continuous variables
%                   (Artelys Knitro) -- not yet implemented
%       ALG ('opf.ac') : (optional) 'opf.ac' or any valid value of
%               OPT.alg for NLPS_MASTER. The first option indicates that
%               it should be taken from MPOPT.opf.ac.solver.
%
%   Output:
%       NLPOPT : an options struct for use by NLPS_MASTER and friends
%
% See also nlps_master, mpoption.

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
    model = 'NLP';
else
    model = upper(model);
end

%% get ALG from mpopt, if necessary
switch alg
    case {'opf.ac', ''}
        alg = upper(mpopt.opf.ac.solver);
    otherwise
        alg = upper(alg);
end

%% default solver
if strcmp(alg, 'DEFAULT')
    alg = 'MIPS';
    if model(1) == 'M'  %% mixed integer
        if have_feature('knitromatlab')
            alg = 'KNITRO';
        else
            error('mpopt2nlpopt: Sorry, no solver available for %s models', model);
        end
    end
end

%% create NLPS_MASTER options struct
nlpopt = struct('alg', alg, 'verbose', mpopt.verbose);
switch alg
    case 'MIPS'
        %% set up options
        nlpopt.mips_opt = mpopt.mips;
        if nlpopt.mips_opt.feastol == 0      %% = MPOPT.opf.violation by default
            nlpopt.mips_opt.feastol = mpopt.opf.violation;
        end
        if ~isfield(nlpopt.mips_opt, 'cost_mult') || isempty(nlpopt.mips_opt.cost_mult)
            nlpopt.mips_opt.cost_mult = 1e-4;
        end
    case 'FMINCON'
        %% basic optimset options needed for fmincon
        nlpopt.fmincon_opt = mpopt.fmincon;
        nlpopt.fmincon_opt.opts.TolCon = mpopt.opf.violation;
    case 'IPOPT'
        nlpopt.ipopt_opt = ipopt_options([], mpopt);
    case 'KNITRO'
        nlpopt.knitro_opt               = mpopt.knitro;
        nlpopt.knitro_opt.opts.feastol  = mpopt.opf.violation;
%         nlpopt.knitro_opt.opts.TolCon   = mpopt.opf.violation;
end
