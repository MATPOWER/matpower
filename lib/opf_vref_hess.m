function d2Vref = opf_vref_hess(x, lam, mpc, refs, mpopt)
%OPF_VREF_HESS  Evaluates Hessian of voltage angle reference.
%   D2VREF = OPF_VREF_HESS(X, LAMBDA, MPC, REFS, MPOPT)
%
%   Hessian evaluation function for voltage angle reference.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPC : MATPOWER case struct
%     REFS : reference vector
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VREF : Hessian of voltage angle reference.
%
%   Example:
%       d2Vref = opf_vref_hess(x, lambda, mpc, refs, mpopt);
%
%   See also OPF_VREF_FCN.

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% unpack data
[Vr, Vi] = deal(x{:});
nb = length(Vr);

%%----- evaluate Hessian of voltage angle constraints -----
VVr = Vr(refs);
VVi = Vi(refs);
VVr2 = VVr.^2;
VVi2 = VVi.^2;
lamVm4 = lam ./ (VVr2 + VVi2).^2;

d2Vref_rr = sparse(refs, refs, 2 * lamVm4 .*  VVr .* VVi,   nb, nb);
d2Vref_ri = sparse(refs, refs,     lamVm4 .* (VVi2 - VVr2), nb, nb);
d2Vref_ir =  d2Vref_ri;
d2Vref_ii = -d2Vref_rr;

%% construct Hessian
d2Vref = [  d2Vref_rr d2Vref_ri;
            d2Vref_ir d2Vref_ii   ];
