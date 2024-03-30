function d2Veq = opf_veq_hess(x, lambda, mpc, idx, mpopt)
% opf_veq_hess - Evaluates Hessian of voltage magnitude equality constraint.
% ::
%
%   D2VEQ = OPF_VEQ_HESS(X, LAMBDA, MPC, IDX, MPOPT)
%
%   Hessian evaluation function for voltage magnitudes.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPC : MATPOWER case struct
%     IDX : index of buses whose voltage magnitudes should be fixed
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VEQ : Hessian of voltage magnitudes.
%
%   Example:
%       d2Veq = opf_veq_hess(x, lambda, mpc, idx, mpopt);
%
% See also opf_veq_fcn.

%   MATPOWER
%   Copyright (c) 2018-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses
n = length(idx);            %% number of buses with fixed voltage magnitudes

%%----- evaluate Hessian of voltage magnitude constraints -----
dlam = sparse(idx, idx, 2*lambda, nb, nb);
zz = sparse(nb, nb);

%% construct Hessian
d2Veq =  [dlam zz; zz dlam];
