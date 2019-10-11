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
%   See https://matpower.org for more info.

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses
n = length(idx);            %% number of buses with voltage limits

%%----- evaluate Hessian of voltage limit constraints -----
nlam = length(lambda) / 2;
if nlam
    lamVmin = lambda(1:nlam);
    lamVmax = lambda((1:nlam)+nlam);
else    %% keep dimensions of empty matrices/vectors compatible
    lamVmin = zeros(0,1);
    lamVmax = zeros(0,1);
end

dlamVmin = sparse(idx, idx, 2*lamVmin, nb, nb);
dlamVmax = sparse(idx, idx, 2*lamVmax, nb, nb);
zz = sparse(nb, nb);

%% construct Hessian
d2Vlims =  -[dlamVmin zz; zz dlamVmin] + ...
            [dlamVmax zz; zz dlamVmax];
