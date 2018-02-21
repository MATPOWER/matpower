function d2G = opf_current_balance_hess(x, lambda, mpc, Ybus, mpopt)
%OPF_CURRENT_BALANCE_HESS  Evaluates Hessian of current balance constraints.
%   D2G = OPF_CURRENT_BALANCE_HESS(X, LAMBDA, OM, YBUS, MPOPT)
%
%   Hessian evaluation function for AC real and imaginary current balance
%   constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of current balance constraints.
%
%   Example:
%       d2G = opf_current_balance_hess(x, lambda, mpc, Ybus, mpopt);
%
%   See also OPF_CURRENT_BALANCE_FCN.

%%----- initialize -----    
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% unpack data
[baseMVA, bus, gen] = deal(mpc.baseMVA, mpc.bus, mpc.gen);

if mpopt.opf.v_cartesian
    [Vi, Vr, Pg, Qg] = deal(x{:});  %% unpack data
    %% reconstruct V
    V = Vr + 1j* Vi;  
    %% problem dimensions
    nb = length(Vi);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections    
else
    [Va, Vm, Pg, Qg] = deal(x{:});  %% unpack data
    %% reconstruct V
    V = Vm .* exp(1j * Va);
    %% problem dimensions
    nb = length(Va);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections    
end
%% put Pg, Qg back in gen
    gen(:,PG) = Pg*baseMVA;
    gen(:,QG) = Qg*baseMVA;
    
%% rebuild Sbus    
    Sbus = makeSbus(baseMVA, bus, gen);
        
%% ----- evaluate Hessian of power balance constraints -----
    nlam = length(lambda) / 2;
    lamP = lambda(1:nlam);
    lamQ = lambda((1:nlam)+nlam);
    Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng);    
    
if mpopt.opf.v_cartesian
    [Grii, Grir, Grri, Grrr] = d2Imis_dV2_C(Sbus, Ybus, V, lamP);
    [Giii, Giir, Giri, Girr] = d2Imis_dV2_C(Sbus, Ybus, V, lamQ);
    [Gr_SV, Gr_VS] = d2Imis_dVdSg_C(V, lamP, Cg);
    [Gi_SV, Gi_VS] = d2Imis_dVdSg_C(V, lamQ, Cg);  
    
    %% construct Hessian
    d2G = [
        real([Grii Grir; Grri Grrr]) + imag([Giii Giir; Giri Girr]), real(Gr_VS) + imag(Gi_VS);
        real(Gr_SV) + imag(Gi_SV), sparse(2*ng, 2*ng)
        ];       
else
    [Graa, Grav, Grva, Grvv] = d2Imis_dV2_P(Sbus, Ybus, V, lamP);
    [Giaa, Giav, Giva, Givv] = d2Imis_dV2_P(Sbus, Ybus, V, lamQ);
    [Gr_SV, Gr_VS] = d2Imis_dVdSg_P(V, lamP, Cg);
    [Gi_SV, Gi_VS] = d2Imis_dVdSg_P(V, lamQ, Cg);

    %% construct Hessian
    d2G = [
        real([Graa Grav; Grva Grvv]) + imag([Giaa Giav; Giva Givv]), real(Gr_VS) + imag(Gi_VS);
        real(Gr_SV) + imag(Gi_SV), sparse(2*ng, 2*ng)
        ];   
end
end


