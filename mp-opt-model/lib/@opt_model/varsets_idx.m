function kk = varsets_idx(om, vs)
% varsets_idx - Returns a vector of indices into X specified by VARSETS
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_variable.varsets_idx instead, as
%       in ``om.var.varsets_idx(...)``.
%
% ::
%
%   K = OM.VARSETS_IDX(VARSETS)
%
%   Returns a vector of indices into X corresponding to the variable
%   sets specified by VARSETS.
%
% See also varsets_x.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

kk = om.var.varsets_idx(vs);
