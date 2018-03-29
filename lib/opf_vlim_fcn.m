function [Vlims, dVlims] = opf_vlim_fcn(x, mpc, mpopt)
%OPF_VLIM_FCN  Evaluates voltage magnitudes and their gradients.
%   [Vlims, dVlims] = OPF_VLIM_FCN(X, mpc, ref, MPOPT)
%
%   Computes the voltage magnitudes using real and imaginary part of complex voltage for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     VLIMS  : vector of voltage magnitudes
%     DVLIMS : (optional) magnitude gradients
%
%   Examples:
%       Vlims = opf_vlim_fcn(x, mpc, mpopt);
%       [Vlims, dVlims] = opf_vlim_fcn(x, mpc, mpopt);
%
%   See also OPF_VLIM_HESS

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Baljinnyam Sereeter, Delft University of Technology
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses

%% compute voltage magnitude
Vm = sqrt(Vr.^2 + Vi.^2);
Vlims = [ mpc.bus(:, VMIN) - Vm;
          Vm - mpc.bus(:, VMAX) ];

if nargout > 1
    %% compute partials of voltage magnitude w.r.t Vr and Vi
    dVm_dVr = sparse(1:nb, 1:nb, Vr./Vm, nb, nb);
    dVm_dVi = sparse(1:nb, 1:nb, Vi./Vm, nb, nb);
    dVlims = [ -dVm_dVr -dVm_dVi;   %% Vlims w.r.t Vr, Vi
                dVm_dVr  dVm_dVi  ];
end