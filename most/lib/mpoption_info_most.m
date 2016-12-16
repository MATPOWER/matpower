function opt = mpoption_info_most(selector)
%MPOPTION_INFO_MOST  Returns MATPOWER option info for MOST.
%
%   DEFAULT_OPTS = MPOPTION_INFO_MOST('D')
%   VALID_OPTS   = MPOPTION_INFO_MOST('V')
%   EXCEPTIONS   = MPOPTION_INFO_MOST('E')
%
%   Returns a structure for MOST options for MATPOWER containing ...
%   (1) default options,
%   (2) valid options, or
%   (3) NESTED_STRUCT_COPY exceptions for setting options
%   ... depending on the value of the input argument.
%
%   This function is used by MPOPTION to set default options, check validity
%   of option names or modify option setting/copying behavior for this
%   subset of optional MATPOWER options.
%
%   See also MPOPTION.

%   MOST
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    selector = 'D';
end
if have_fcn('most')
    switch upper(selector)
        case {'D', 'V'}     %% default and valid options
            opt = struct(...
                'most',   struct(...
                    'build_model',                  1, ...          %% was md.CreateQP
                    'solve_model',                  1, ...          %% was md.Solve
                    'resolve_new_cost',             0, ...          %% was md.ReSolveNewCoordCost
                    'dc_model',                     1, ...          %% was md.DCMODEL = []
                    'fixed_res',                    -1, ...         %% was md.IncludeFixedReserves = []
                    'q_coordination',               0, ...          %% was md.QCoordination
                    'security_constraints',         -1, ...         %% was md.SecurityConstrained = []
                    'storage',                      struct(...
                        'terminal_target',              -1, ...     %% was md.Storage.ForceExpectedTerminalStorage
                        'cyclic',                       0), ...     %% was md.Storage.ForceCyclicStorage
                    'uc',                           struct(...
                        'run',                          -1, ...     %% was missing
                        'cyclic',                       0), ...     %% was md.UC.CyclicCommitment
                    'alpha',                        0, ...          %% was md.alpha
                    'solver',                       'DEFAULT', ...  %% was md.QP.opt.alg
                    'skip_prices',                  0, ...          %% was md.QP.opt.skip_prices
                    'price_stage_warn_tol',         1e-7 ...
                ) ...
            );
        case 'E'            %% exceptions used by nested_struct_copy() for applying
            opt = struct([]);   %% no exceptions
        otherwise
            error('mpoption_info_most: ''%s'' is not a valid input argument', selector);
    end
else
    opt = struct([]);       %% MOST is not available
end
