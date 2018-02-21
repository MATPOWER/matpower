function [g, dg] = opf_current_balance_fcn(x, mpc, Ybus, mpopt)
%OPF_CURRENT_BALANCE_FCN  Evaluates AC current balance constraints and their gradients.
%   [G, DG] = OPF_CURRENT_BALANCE_FCN(X, OM, YBUS, MPOPT)
%
%   Computes the real or imaginary current balance equality constraints for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (real/imaginary current balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_current_balance_fcn(x, mpc, Ybus, mpopt);
%       [g, dg] = opf_current_balance_fcn(x, mpc, Ybus, mpopt);
%
%   See also OPF_CURRENT_BALANCE_HESS

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% unpack data
[baseMVA, bus, gen] = deal(mpc.baseMVA, mpc.bus, mpc.gen);

if mpopt.opf.v_cartesian
    [Vi, Vr, Pg, Qg] = deal(x{:});
    %% reconstruct V
    V = Vr + 1j* Vi;  
    %% problem dimensions
    nb = length(Vi);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections        
else
    [Va, Vm, Pg, Qg] = deal(x{:});
    %% reconstruct V
    V = Vm .* exp(1j * Va);
    %% problem dimensions
    nb = length(Va);            %% number of buses
    ng = length(Pg);            %% number of dispatchable injections    
end 

%% ----- evaluate constraints -----
%% put Pg, Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

%% rebuild Sbus
Sbus = makeSbus(baseMVA, bus, gen);  %% net injected power in p.u.

%% evaluate complex current balance mismatches
mis = Ybus*V - conj(Sbus./V);

%% assemble active and reactive power balance constraints
g = [ real(mis);    %% real current mismatch
      imag(mis) ];  %% imaginary current mismatch

%%----- evaluate constraint gradients -----
if nargout > 1
    %% compute partials of injected bus powers
    Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng);
    InvConjV = sparse(1:nb, 1:nb, 1./conj(V), nb, nb);
    dImis_dPg = -InvConjV*Cg;       % dImis w.r.t. Pg
    dImis_dQg = 1j*InvConjV*Cg;     % dImis w.r.t. Qg

    if mpopt.opf.v_cartesian
        [dImis_dVr, dImis_dVi] = dImis_dV_C(Sbus, Ybus, V);          % w.r.t. V
        dg = [
            real([dImis_dVi dImis_dVr dImis_dPg dImis_dQg]);   %% Ir mismatch w.r.t Vi, Vr, Pg, Qg
            imag([dImis_dVi dImis_dVr dImis_dPg dImis_dQg]) ;  %% Ii mismatch w.r.t Vi, Vr, Pg, Qg
        ];        
    else
        [dImis_dVm, dImis_dVa] = dImis_dV_P(Sbus, Ybus, V);          % w.r.t. V
        dg = [
            real([dImis_dVa dImis_dVm dImis_dPg dImis_dQg]);  %% Ir mismatch w.r.t Va, Vm, Pg, Qg
            imag([dImis_dVa dImis_dVm dImis_dPg dImis_dQg]);  %% Ii mismatch w.r.t Va, Vm, Pg, Qg
        ];
    end
end
