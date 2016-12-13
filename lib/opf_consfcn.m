function [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il, varargin)
%OPF_CONSFCN  Evaluates nonlinear constraints and their Jacobian for OPF.
%   [H, G, DH, DG] = OPF_CONSFCN(X, OM, YBUS, YF, YT, MPOPT, IL)
%
%   Constraint evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON. Computes constraint vectors and their
%   gradients.
%
%   Inputs:
%     X : optimization vector
%     OM : OPF model object
%     YBUS : bus admittance matrix
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     MPOPT : MATPOWER options struct
%     IL : (optional) vector of branch indices corresponding to
%          branches with flow limits (all others are assumed to be
%          unconstrained). The default is [1:nl] (all branches).
%          YF and YT contain only the rows corresponding to IL.
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          limit^2 - flow^2, where the flow can be apparent power
%          real power or current, depending on value of
%          opf.flow_lim in MPOPT (only for constrained lines)
%     G  : vector of equality constraint values (power balances)
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       [h, g] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt);
%       [h, g, dh, dg] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il);
%
%   See also OPF_COSTFCN, OPF_HESSFCN.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   Ray Zimmerman, PSERC Cornell and Shrirang Abhyankar
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
vv = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types

%% set default constrained lines
if nargin < 7
    il = (1:nl);            %% all lines have limits by default
end
nl2 = length(il);           %% number of constrained lines

%% grab Pg & Qg
Pg = x(vv.i1.Pg:vv.iN.Pg);  %% active generation in p.u.
Qg = x(vv.i1.Qg:vv.iN.Qg);  %% reactive generation in p.u.

%% put Pg & Qg back in gen
gen(:, PG) = Pg * baseMVA;  %% active generation in MW
gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr
 
%% ----- evaluate constraints -----
%% reconstruct V
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
V = Vm .* exp(1j * Va);

%% rebuild Sbus
Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.

%% evaluate power flow equations
mis = V .* conj(Ybus * V) - Sbus;

%%----- evaluate constraint function values -----
%% first, the equality constraints (power flow)
g = [ real(mis);            %% active power mismatch for all buses
      imag(mis) ];          %% reactive power mismatch for all buses

%% then, the inequality constraints (branch flow limits)
if nl2 > 0
  flow_max = (branch(il, RATE_A)/baseMVA).^2;
  flow_max(flow_max == 0) = Inf;
  if upper(mpopt.opf.flow_lim(1)) == 'I'    %% current magnitude limit, |I|
    If = Yf * V;
    It = Yt * V;
    h = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
          It .* conj(It) - flow_max ];  %% branch current limits (to bus)
  else
    %% compute branch power flows
    Sf = V(branch(il, F_BUS)) .* conj(Yf * V);  %% complex power injected at "from" bus (p.u.)
    St = V(branch(il, T_BUS)) .* conj(Yt * V);  %% complex power injected at "to" bus (p.u.)
    if upper(mpopt.opf.flow_lim(1)) == 'P'  %% active power limit, P (Pan Wei)
      h = [ real(Sf).^2 - flow_max;         %% branch real power limits (from bus)
            real(St).^2 - flow_max ];       %% branch real power limits (to bus)
    else                                    %% apparent power limit, |S|
      h = [ Sf .* conj(Sf) - flow_max;      %% branch apparent power limits (from bus)
            St .* conj(St) - flow_max ];    %% branch apparent power limits (to bus)
    end
  end
else
  h = zeros(0,1);
end

%%----- evaluate partials of constraints -----
if nargout > 2
  %% index ranges
  iVa = vv.i1.Va:vv.iN.Va;
  iVm = vv.i1.Vm:vv.iN.Vm;
  iPg = vv.i1.Pg:vv.iN.Pg;
  iQg = vv.i1.Qg:vv.iN.Qg;

  %% compute partials of injected bus powers
  [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);           %% w.r.t. V
  [dummy, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm);
  dSbus_dVm = dSbus_dVm - neg_dSd_dVm;
  neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
                                                        %% Qbus w.r.t. Qg
  
  %% construct Jacobian of equality (power flow) constraints and transpose it
  dg = sparse(2*nb, nxyz);
  dg(:, [iVa iVm iPg iQg]) = [
    real([dSbus_dVa dSbus_dVm]) neg_Cg sparse(nb, ng);  %% P mismatch w.r.t Va, Vm, Pg, Qg
    imag([dSbus_dVa dSbus_dVm]) sparse(nb, ng) neg_Cg;  %% Q mismatch w.r.t Va, Vm, Pg, Qg
  ];
  dg = dg';

  if nl2 > 0
    %% compute partials of Flows w.r.t. V
    if upper(mpopt.opf.flow_lim(1)) == 'I'  %% current
      [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V);
    else                            %% power
      [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, V);
    end
    if upper(mpopt.opf.flow_lim(1)) == 'P'  %% real part of flow (active power)
      dFf_dVa = real(dFf_dVa);
      dFf_dVm = real(dFf_dVm);
      dFt_dVa = real(dFt_dVa);
      dFt_dVm = real(dFt_dVm);
      Ff = real(Ff);
      Ft = real(Ft);
    end
  
    %% squared magnitude of flow (of complex power or current, or real power)
    [df_dVa, df_dVm, dt_dVa, dt_dVm] = ...
            dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
  
    %% construct Jacobian of inequality (branch flow) constraints & transpose
    dh = sparse(2*nl2, nxyz);
    dh(:, [iVa iVm]) = [
      df_dVa, df_dVm;                     %% "from" flow limit
      dt_dVa, dt_dVm;                     %% "to" flow limit
    ];
    dh = dh';
  else
    dh = sparse(nxyz, 0);
  end
end
