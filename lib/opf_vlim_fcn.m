function [Vm, dVm] = opf_vlim_fcn(x,mpc, ref, mpopt)

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