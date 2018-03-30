function [VaDif, dVaDif] = opf_branch_ang_fcn(x, Aang, lang, uang, iang, mpopt);
%OPF_BRANCH_ANG_FCN  Evaluates branch angle difference constraints and gradients.
%   [VADIF, DVADIF] = OPF_BRANCH_ANG_FCN(X, AANG, LANG, UANG, IANG, MPOPT);
%
%   Computes the lower and upper constraints on branch angle differences
%   for voltages in cartesian coordinates. Computes constraint vectors and
%   their gradients. The constraints are of the form:
%       Aang * Va >= lang
%       Aang * Va <= uang
%   where Va is the voltage angle, a non-linear function of the Vr and Vi.
%
%   Inputs:
%     X : optimization vector
%     AANG : constraint matrix, see MAKEAANG
%     LANG : lower bound vector, see MAKEAANG
%     UANG : upper bound vector, see MAKEAANG
%     IANG : index vector of branches corresponding to rows of AANG, LANG, UANG
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     VADIF  : constraint vector [ lang - Aang * Va; Aang * Va - uang ]
%     DVADIF : (optional) constraint gradients
%
%   Examples:
%       VaDif = opf_branch_ang_fcn(x, Aang, lang, uang, iang, mpopt);
%       [VaDif, dVaDif] = opf_branch_ang_fcn(x, Aang, lang, uang, iang, mpopt);
%
%   See also OPF_BRANCH_ANG_HESS

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

%% problem dimensions
nb = length(Vi);            %% number of buses

%% compute branch angle difference
Va = angle(Vr + 1j* Vi);
Ax = Aang * Va;
VaDif = [ lang - Ax;
          Ax - uang ];

if nargout > 1
    %% compute partials of branch angle difference w.r.t Vr and Vi
    Vm2 = Vr.^2 + Vi.^2;
    AangdVa_dVr = Aang * sparse(1:nb, 1:nb, -Vi./Vm2, nb, nb);
    AangdVa_dVi = Aang * sparse(1:nb, 1:nb,  Vr./Vm2, nb, nb);
    dVaDif = [ -AangdVa_dVr -AangdVa_dVi;   %% VaDif w.r.t Vr, Vi
                AangdVa_dVr  AangdVa_dVi ];
end
