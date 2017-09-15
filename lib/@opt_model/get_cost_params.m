function cp = get_cost_params(om, name, idx)
%GET_COST_PARAMS  Returns the cost parameter struct for user-defined costs.
%
%   -----  DEPRECATED - use PARAMS_LEGACY_COST instead  -----
%
%   CP = OM.GET_COST_PARAMS()
%   CP = OM.GET_COST_PARAMS(NAME)
%   CP = OM.GET_COST_PARAMS(NAME, IDX)
%
%   Returns the full cost parameter struct for all user-defined
%   costs that incorporates all of the named cost sets added via ADD_COSTS,
%   or, if a NAME is provided it returns the cost struct corresponding to
%   the named set of cost rows (N still has full number of columns).
%
%   The cost parameters are returned in a struct with the following fields:
%       N      - nw x nx sparse matrix
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   See also OPT_MODEL, ADD_COSTS, COMPUTE_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

cp = om.params_legacy_cost();

if nargin > 1
    if getN(om, 'cost', name)
        if nargin < 3 || isempty(idx)
            if prod(size(om.cost.idx.i1.(name))) == 1
                i1 = om.cost.idx.i1.(name);
                iN = om.cost.idx.iN.(name);
            else
                error('@opt_model/get_cost_params: cost set ''%s'' requires an idx arg', name);
            end
        else
            %% calls to substruct() are relatively expensive, so we pre-build the
            %% structs for addressing cell and numeric array fields
            %% sn = substruct('.', name, '()', idx);
            sn = struct('type', {'.', '()'}, 'subs', {name, idx});  %% num array field
            i1 = subsref(om.cost.idx.i1, sn);
            iN = subsref(om.cost.idx.iN, sn);
        end
        cp.N  = cp.N(i1:iN,:);
        cp.Cw = cp.Cw(i1:iN);
        cp.H  = cp.H(i1:iN,i1:iN);
        cp.dd = cp.dd(i1:iN);
        cp.rh = cp.rh(i1:iN);
        cp.kk = cp.kk(i1:iN);
        cp.mm = cp.mm(i1:iN);
    end
end
