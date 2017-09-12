function [A, l, u] = linear_constraints(om)
%LINEAR_CONSTRAINTS  Builds and returns the full set of linear constraints.
%
%   -----  DEPRECATED - Please use PARAMS_LIN_CONSTRAINT instead  -----
%
%   [A, L, U] = OM.LINEAR_CONSTRAINTS()
%   Builds the full set of linear constraints based on those added by
%   ADD_LIN_CONSTRAINT.
%
%       L <= A * x <= U
%
%   Example:
%       [A, l, u] = om.linear_constraints();
%
%   See also OPT_MODEL, ADD_LIN_CONSTRAINT, PARAMS_LIN_CONSTRAINT.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[A, l, u] = om.params_lin_constraint();
