function d2Vlims = opf_vlim_hess(x, lambda, mpc, mpopt)
%OPF_VLIM_HESS  Evaluates Hessian of voltage magnitudes.
%   D2VLIMS = OPF_VLIM_HESS(X, LAMBDA, MPC, MPOPT)
%
%   Hessian evaluation function for voltage magnitudes.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPC : MATPOWER case struct
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VLIMS : Hessian of voltage magnitudes.
%
%   Example:
%       d2Vlims = opf_vlim_hess(x, lambda, mpc, mpopt);
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

%% compute voltage magnitude cubed
Vr2 = Vr.^2;
Vi2 = Vi.^2;
VrVi = Vr .* Vi;
Vm3 = (Vr2 + Vi2).^(3/2);   %% Vm.^3;

%%----- evaluate Hessian of voltage limit constraints -----
nlam = length(lambda) / 2;
if nlam
    lamVmax = lambda(1:nlam);
    lamVmin = lambda((1:nlam)+nlam);
else    %% keep dimensions of empty matrices/vectors compatible
    lamVmax = zeros(0,1);
    lamVmin = zeros(0,1);
end

lamVmax_over_Vm3 = lamVmax ./ Vm3;
lamVmin_over_Vm3 = lamVmin ./ Vm3;

Vmax_rr = sparse(1:nb, 1:nb,  Vi2  .* lamVmax_over_Vm3, nb, nb);
Vmax_ri = sparse(1:nb, 1:nb, -VrVi .* lamVmax_over_Vm3, nb, nb);
Vmax_ir = Vmax_ri;
Vmax_ii = sparse(1:nb, 1:nb,  Vr2  .* lamVmax_over_Vm3, nb, nb);

Vmin_rr = sparse(1:nb, 1:nb,  Vi2  .* lamVmin_over_Vm3, nb, nb);
Vmin_ri = sparse(1:nb, 1:nb, -VrVi .* lamVmin_over_Vm3, nb, nb);
Vmin_ir = Vmin_ri;
Vmin_ii = sparse(1:nb, 1:nb,  Vr2  .* lamVmin_over_Vm3, nb, nb);

%% construct Hessian
d2Vlims =  -[Vmin_rr Vmin_ri; Vmin_ir Vmin_ii] + ...
            [Vmax_rr Vmax_ri; Vmax_ir Vmax_ii];
