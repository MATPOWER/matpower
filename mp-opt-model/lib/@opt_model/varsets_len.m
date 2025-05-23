function nv = varsets_len(om, vs)
% varsets_len - Returns the total number of variables in VARSETS
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_variable.varsets_len instead, as
%       in ``om.var.varsets_len(...)``.
%
% ::
%
%   NV = OM.VARSETS_LEN(VARSETS)
%
%   Returns the total number of elements in the optimization sub-vector
%   specified by VARSETS.
%
% See also varsets_cell2struct.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

nv = om.var.varsets_len(vs);
