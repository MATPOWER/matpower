function d2Vm = opf_vlim_hess(x, lam, ref, mpopt)
%OPF_VLIM_HESS  Evaluates Hessian of voltage magnitudes.
%   D2VM = OPF_VLIM_HESS(X, LAMBDA, REF, MPOPT)
%
%   Hessian evaluation function for voltage magnitudes.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2VM : Hessian of voltage magnitudes.
%
%   Example:
%       d2VM = opf_vlim_hess(x, lambda, ref, mpopt);
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

% nlam = length(lambda) / 2;
% lamP = lambda(1:nlam);
% lamQ = lambda((1:nlam)+nlam);

%% ----- evaluate constraint gradients -----
diagVi       = sparse(1:n, 1:n, Vi, n, n);    
diagVr       = sparse(1:n, 1:n, Vr, n, n);    
diagVm1lam   = sparse(1:n, 1:n, lam./Vm, n, n);    
diagVilamVm3 = sparse(1:n, 1:n, (Vi.*lam)./(Vm.^3), n, n);    
diagVrlamVm3 = sparse(1:n, 1:n, (Vr.*lam)./(Vm.^3), n, n);    

Vm_ii = diagVm1lam - diagVi*diagVilamVm3;
Vm_ir = - diagVi*diagVrlamVm3 ;
Vm_ri = Vm_ir;
Vm_rr = diagVm1lam - diagVr*diagVrlamVm3;
%% construct Hessian
d2Vm = [ Vm_ii Vm_ir; Vm_ri Vm_rr];    

end