function xx = varsets_x(om, x, varargin)
% varsets_x - Returns a cell array of sub-vectors of X specified by VARSETS
%
% .. note::
%    .. deprecated:: 5.0 Please use mp.sm_variable.varsets_x instead, as
%       in ``om.var.varsets_x(...)``.
%
% ::
%
%   X = OM.VARSETS_X(X, VARSETS)
%   X = OM.VARSETS_X(X, VARSETS, 'vector')
%
%   Returns a cell array of sub-vectors of X specified by VARSETS, or
%   the full optimization vector X, if VARSETS is empty.
%
%   If a 3rd argument is present (value is ignored) the returned value is
%   a single numeric vector with the individual components stacked vertically.
%
% See also varsets_len.

%   MP-Opt-Model
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

xx = om.var.varsets_x(x, varargin{:});
