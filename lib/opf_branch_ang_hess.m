function d2VaDif = opf_branch_ang_hess(x, lambda, Aang, lang, uang, iang, mpopt)
%OPF_BRANCH_ANG_HESS  Evaluates Hessian of branch angle difference constraints.
%   D2VADIF = OPF_BRANCH_ANG_HESS(X, LAMBDA, AANG, LANG, UANG, IANG, MPOPT)
%
%   Hessian evaluation function for branch angle difference constraints
%   for voltages in cartesian coordinates.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on branch angle
%           difference constraints, lower, then upper
%     AANG : constraint matrix, see MAKEAANG
%     LANG : lower bound vector, see MAKEAANG
%     UANG : upper bound vector, see MAKEAANG
%     IANG : index vector of branches corresponding to rows of AANG, LANG, UANG
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VADIF : Hessian of branch angle difference constraints.
%
%   Example:
%       d2VaDif = opf_branch_ang_hess(x, lambda, Aang, lang, uang, iang, mpopt);
%
%   See also OPF_BRANCH_ANG_FCN.

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% unpack data
[Vr, Vi] = deal(x{:});
nb = length(Vr);

%%----- evaluate Hessian of branch angle difference constraints -----
nlam = length(lambda) / 2;
if nlam
    lamL = lambda(1:nlam);
    lamU = lambda((1:nlam)+nlam);
else
    lamL = zeros(0,1);
    lamU = zeros(0,1);
end

Vr2 = Vr.^2;
Vi2 = Vi.^2;

lamL_Vm4 = (Aang' * lamL) ./ (Vr2 + Vi2).^2;
lamU_Vm4 = (Aang' * lamU) ./ (Vr2 + Vi2).^2;

VaDifL_rr = sparse(1:nb, 1:nb, 2 * lamL_Vm4 .*  Vr .* Vi,   nb, nb);
VaDifL_ri = sparse(1:nb, 1:nb,     lamL_Vm4 .* (Vi2 - Vr2), nb, nb);
VaDifL_ir =  VaDifL_ri;
VaDifL_ii = -VaDifL_rr;

VaDifU_rr = sparse(1:nb, 1:nb, 2 * lamU_Vm4 .*  Vr .* Vi,   nb, nb);
VaDifU_ri = sparse(1:nb, 1:nb,     lamU_Vm4 .* (Vi2 - Vr2), nb, nb);
VaDifU_ir =  VaDifU_ri;
VaDifU_ii = -VaDifU_rr;

%% construct Hessian
d2VaDif = -[ VaDifL_rr VaDifL_ri;
             VaDifL_ir VaDifL_ii ] + ...
           [ VaDifU_rr VaDifU_ri;
             VaDifU_ir VaDifU_ii ];
