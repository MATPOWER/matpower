function d2Vlims = opf_vlim_hess(x, lambda, mpc, idx, mpopt)
% opf_vlim_hess - Evaluates Hessian of voltage magnitudes.
% ::
%
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
% See also opf_vlim_fcn.

%   MATPOWER
%   Copyright (c) 2018-2024, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses
n = length(idx);            %% number of buses with voltage limits

%%----- evaluate Hessian of voltage limit constraints -----
nlam = length(lambda) / 2;
if nlam
    lam = lambda(nlam+1:2*nlam) - lambda(1:nlam);   %% lamVmax - lamVmin
else    %% keep dimensions of empty matrices/vectors compatible
    lam = zeros(0,1);
end

dlam = sparse(idx, idx, 2*lam, nb, nb);
zz = sparse(nb, nb);

%% construct Hessian
d2Vlims = [dlam zz; zz dlam];
