function d2Vlims = opf_vlim_hess(x, lambda, mpc, idx, mpopt)
%OPF_VLIM_HESS  Evaluates Hessian of voltage magnitudes.
%   D2VLIMS = OPF_VLIM_HESS(X, LAMBDA, MPC, IDX, MPOPT)
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
%     D2VLIMS : Hessian of voltage magnitudes.
%
%   Example:
%       d2Vlims = opf_vlim_hess(x, lambda, mpc, idx, mpopt);
%
%   See also OPF_VLIM_FCN.

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
n = length(idx);            %% number of buses with voltage limits

%% compute voltage magnitude cubed
Vr2 = Vr(idx).^2;
Vi2 = Vi(idx).^2;
VrVi = Vr(idx) .* Vi(idx);
Vm3 = (Vr2 + Vi2).^(3/2);   %% Vm.^3;

%%----- evaluate Hessian of voltage limit constraints -----
nlam = length(lambda) / 2;
if nlam
    lamVmin = lambda(1:nlam);
    lamVmax = lambda((1:nlam)+nlam);
else    %% keep dimensions of empty matrices/vectors compatible
    lamVmin = zeros(0,1);
    lamVmax = zeros(0,1);
end

lamVmin_over_Vm3 = lamVmin ./ Vm3;
lamVmax_over_Vm3 = lamVmax ./ Vm3;

Vmin_rr = sparse(idx, idx,  Vi2  .* lamVmin_over_Vm3, nb, nb);
Vmin_ri = sparse(idx, idx, -VrVi .* lamVmin_over_Vm3, nb, nb);
Vmin_ir = Vmin_ri;
Vmin_ii = sparse(idx, idx,  Vr2  .* lamVmin_over_Vm3, nb, nb);

Vmax_rr = sparse(idx, idx,  Vi2  .* lamVmax_over_Vm3, nb, nb);
Vmax_ri = sparse(idx, idx, -VrVi .* lamVmax_over_Vm3, nb, nb);
Vmax_ir = Vmax_ri;
Vmax_ii = sparse(idx, idx,  Vr2  .* lamVmax_over_Vm3, nb, nb);

%% construct Hessian
d2Vlims =  -[Vmin_rr Vmin_ri; Vmin_ir Vmin_ii] + ...
            [Vmax_rr Vmax_ri; Vmax_ir Vmax_ii];
