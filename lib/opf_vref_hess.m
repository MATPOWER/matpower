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
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%% unpack data
[Vr, Vi] = deal(x{:});

%% ----- evaluate constraint gradients -----
diagVrVm4   = diag(Vr./(Vr(refs).^4 + 2*(Vr(refs).^2).*(Vi(refs).^2)+ Vi(refs).^4));
diagViVm4   = diag(Vi./(Vr(refs).^4 + 2*(Vr(refs).^2).*(Vi(refs).^2)+ Vi(refs).^4));
diaglamV2 = diag(lam./(Vr(refs).^2 + Vi(refs).^2));
diagVrlam = diag(2*Vr(refs).*lam);
diagVilam = diag(2*Vi(refs).*lam);

Vref_ii = -diagVrlam*diagViVm4; 
Vref_ir = diaglamV2 - diagVrlam*diagVrVm4;
Vref_ri = -Vref_ir;
Vref_rr = diagVilam*diagVrVm4;
%% construct Hessian
d2Vref = [ Vref_ii Vref_ir; Vref_ri Vref_rr];
end