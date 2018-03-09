function [Vm, dVm] = opf_vlim_fcn(x,mpc, ref, mpopt)
%OPF_VLIM_FCN  Evaluates voltage magnitudes and their gradients.
%   [Vm, dVm] = OPF_VLIM_FCN(X, mpc, ref, MPOPT)
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
%     VM  : vector of voltage magnitudes
%     DVM : (optional) magnitude gradients
%
%   Examples:
%       Vm = opf_vlim_fcn(x, mpc, ref, mpopt);
%       [Vm, dVm] = opf_vlim_fcn(x, mpc, ref, mpopt);
%
%   See also OPF_VLIM_HESS

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses

%% compute voltage magnitude
Vm = sqrt(Vr.^2 + Vi.^2);

%%----- evaluate constraint gradients -----
if nargout > 1
    %% compute partials of voltage magnitue w.r.t Vr and Vi
    dVm_dVr = sparse(1:nb, 1:nb, Vr./Vm, nb, nb);
    dVm_dVi = sparse(1:nb, 1:nb, Vi./Vm, nb, nb);
    dVm = [dVm_dVi dVm_dVr];        %% Vm w.r.t Vi, Vr
end