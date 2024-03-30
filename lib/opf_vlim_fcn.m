function [Vlims, dVlims] = opf_vlim_fcn(x, mpc, idx, mpopt)
% opf_vlim_fcn - Evaluates voltage magnitudes and their gradients.
% ::
%
%   [Vlims, dVlims] = OPF_VLIM_FCN(X, MPC, IDX, MPOPT)
%
%   Computes the voltage magnitudes using real and imaginary part of complex voltage for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     IDX : index of buses whose voltage magnitudes should be fixed
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     VLIMS  : vector of voltage magnitudes
%     DVLIMS : (optional) magnitude gradients
%
%   Examples:
%       Vlims = opf_vlim_fcn(x, mpc, mpopt);
%       [Vlims, dVlims] = opf_vlim_fcn(x, mpc, idx, mpopt);
%
% See also opf_vlim_hess.

%   MATPOWER
%   Copyright (c) 2018-2024, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses
n = length(idx);            %% number of buses with voltage limits

%% compute voltage magnitude
Vm2 = Vr(idx).^2 + Vi(idx).^2;
Vlims = [ mpc.bus(idx, VMIN).^2 - Vm2;
          Vm2 - mpc.bus(idx, VMAX).^2 ];

if nargout > 1
    %% compute partials of voltage magnitude w.r.t Vr and Vi
    dVm_dVr = sparse(1:n, idx, 2 * Vr(idx), n, nb);
    dVm_dVi = sparse(1:n, idx, 2 * Vi(idx), n, nb);
    dVlims = [ -dVm_dVr -dVm_dVi;   %% Vlims w.r.t Vr, Vi
                dVm_dVr  dVm_dVi  ];
end
