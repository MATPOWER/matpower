function [Veq, dVeq] = opf_veq_fcn(x, mpc, idx, mpopt)
%OPF_VEQ_FCN  Evaluates voltage magnitude equality constraint and gradients.
%   [Veq, dVeq] = OPF_VEQ_FCN(X, MPC, IDX, MPOPT)
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
%     VEQ  : vector of voltage magnitudes
%     DVEQ : (optional) magnitude gradients
%
%   Examples:
%       Veq = opf_veq_fcn(x, mpc, mpopt);
%       [Veq, dVeq] = opf_veq_fcn(x, mpc, idx, mpopt);
%
%   See also OPF_VEQ_HESS

%   MATPOWER
%   Copyright (c) 2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
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
n = length(idx);            %% number of buses with fixed voltage magnitudes

%% compute voltage magnitude
Vm = sqrt(Vr(idx).^2 + Vi(idx).^2);
Veq = Vm - mpc.bus(idx, VMAX);

if nargout > 1
    %% compute partials of voltage magnitude w.r.t Vr and Vi
    dVm_dVr = sparse(1:n, idx, Vr(idx)./Vm, n, nb);
    dVm_dVi = sparse(1:n, idx, Vi(idx)./Vm, n, nb);
    dVeq = [ dVm_dVr dVm_dVi ];     %% Vm w.r.t Vr, Vi
end