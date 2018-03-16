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
      
%% ----- evaluate constraint gradients -----
nlam = length(lambda) / 2;
if nlam
    lamUp = lambda(1:nlam);
    lamLow = lambda((1:nlam)+nlam);
else
    lamUp = zeros(0,1);   
    lamLow = zeros(0,1);  
end

%% ----- evaluate constraint gradients -----
diagVi       = sparse(1:n, 1:n, Vi, n, n);    
diagVr       = sparse(1:n, 1:n, Vr, n, n);    
%% for upper limit
    diagVm1lam   = sparse(1:n, 1:n, lamUp./Vm, n, n);    
    diagVilamVm3 = sparse(1:n, 1:n, (Vi.*lamUp)./(Vm.^3), n, n);    
    diagVrlamVm3 = sparse(1:n, 1:n, (Vr.*lamUp)./(Vm.^3), n, n);    

    VlimU_ii = diagVm1lam - diagVi*diagVilamVm3;
    VlimU_ir = - diagVi*diagVrlamVm3 ;
    VlimU_ri = Vm_ir;
    VlimU_rr = diagVm1lam - diagVr*diagVrlamVm3;
%% for lower limit
    diagVm1lam   = sparse(1:n, 1:n, lamLow./Vm, n, n);    
    diagVilamVm3 = sparse(1:n, 1:n, (Vi.*lamLow)./(Vm.^3), n, n);    
    diagVrlamVm3 = sparse(1:n, 1:n, (Vr.*lamLow)./(Vm.^3), n, n);    

    VlimL_ii = diagVm1lam - diagVi*diagVilamVm3;
    VlimL_ir = - diagVi*diagVrlamVm3 ;
    VlimL_ri = Vm_ir;
    VlimL_rr = diagVm1lam - diagVr*diagVrlamVm3;    
%% construct Hessian
d2Vlims = -[ VlimL_ii VlimL_ir; VlimL_ri VlimL_rr] + [ VlimU_ii VlimU_ir; VlimU_ri VlimU_rr];     

end