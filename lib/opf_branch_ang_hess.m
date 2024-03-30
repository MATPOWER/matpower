function d2VaDif = opf_branch_ang_hess(x, lambda, Aang, lang, uang)
% opf_branch_ang_hess - Evaluates Hessian of branch angle difference constraints.
% ::
%
%   D2VADIF = OPF_BRANCH_ANG_HESS(X, LAMBDA, AANG, LANG, UANG)
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
%
%   Outputs:
%     D2VADIF : Hessian of branch angle difference constraints.
%
%   Example:
%       d2VaDif = opf_branch_ang_hess(x, lambda, Aang, lang, uang);
%
% See also opf_branch_ang_fcn.

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
nb = length(Vr);

%%----- evaluate Hessian of branch angle difference constraints -----
nlam = length(lambda) / 2;
if nlam
    lam = lambda(nlam+1:2*nlam) - lambda(1:nlam);   %% lamU - lamL
else
    lam = zeros(0,1);
end

Vr2 = Vr.^2;
Vi2 = Vi.^2;

lam_Vm4 = (Aang' * lam) ./ (Vr2 + Vi2).^2;

VaDif_rr = sparse(1:nb, 1:nb, 2 * lam_Vm4 .*  Vr .* Vi,   nb, nb);
VaDif_ri = sparse(1:nb, 1:nb,     lam_Vm4 .* (Vi2 - Vr2), nb, nb);
VaDif_ir =  VaDif_ri;
VaDif_ii = -VaDif_rr;

%% construct Hessian
d2VaDif = [ VaDif_rr VaDif_ri;
            VaDif_ir VaDif_ii ];
