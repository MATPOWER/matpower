function [VaDif, dVaDif] = opf_branch_ang_fcn(x, Aang, lang, uang, iang, mpopt);

%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses

%% compute branch angle difference
Va = angle(Vr + 1j* Vi);
Ax = Aang * Va;
VaDif = [ lang - Ax;
             Ax - uang ];

if nargout > 1
    %% compute partials of branch angle difference w.r.t Vr and Vi
    Vm2 = Vr.^2 + Vi.^2;
    AangdVa_dVr = Aang*sparse(1:nb, 1:nb, -Vi./Vm2, nb, nb);
    AangdVa_dVi = Aang*sparse(1:nb, 1:nb,  Vr./Vm2, nb, nb);
    dVaDif = [ -[AangdVa_dVr AangdVa_dVi];                      %% VaDif w.r.t Vr, Vi
                [AangdVa_dVr AangdVa_dVi]];
end
