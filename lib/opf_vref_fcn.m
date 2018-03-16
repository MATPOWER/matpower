function [Vref, dVref] = opf_vref_fcn(x, mpc, refs, mpopt)
%OPF_VREF_FCN  Evaluates voltage angle reference and their gradients.
%   [Vref, dVref] = OPF_VREF_FCN(X, mpc, ref, MPOPT)
%
%   Computes the voltage angle reference using real and imaginary part of complex voltage for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     REFS : reference vector
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     VREF  : vector of voltage angle reference
%     DVREF : (optional) angle reference gradients
%
%   Examples:
%       Vref = opf_vref_fcn(x, mpc, refs, mpopt);
%       [Vref, dVref] = opf_vref_fcn(x, mpc, refs, mpopt);
%
%   See also OPF_VREF_HESS

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
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

%% compute voltage angle reference
    Vref = angle(Vr(refs) + 1j* Vi(refs)) - mpc.bus(refs, VA);
    
if nargout > 1
    %% compute partials of voltage angle reference w.r.t Vr and Vi
    Vm2 = Vr(refs).^2 + Vi(refs).^2;
    dVa_dVr = diag(-Vi(refs)./Vm2);
    dVa_dVi = diag(Vr(refs)./Vm2);
    dVref = [dVa_dVr dVa_dVi];                           %% Vref w.r.t Vr, Vi
end    
end