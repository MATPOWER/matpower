function [h, dh] = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt)
%OPF_BRANCH_FLOW_FCN  Evaluates AC branch flow constraints and Jacobian.
%   [H, DH] = OPF_BRANCH_FLOW_FCN(X, OM, YF, YT, IL, MPOPT)
%
%   Active power balance equality constraints for AC optimal power flow.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     IL : vector of branch indices corresponding to branches with
%          flow limits (all others are assumed to be unconstrained).
%          YF and YT contain only the rows corresponding to IL.
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     H  : vector of inequality constraint values (flow limits)
%          where the flow can be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (limit^2 - flow^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow).
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%
%   Examples:
%       h = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%       [h, dh] = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_HESS.

%   MATPOWER
%   Copyright (c) 1996-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% unpack data
lim_type = upper(mpopt.opf.flow_lim(1));
branch = mpc.branch;
if mpopt.opf.v_cartesian
    [Vr, Vi] = deal(x{:});
    V = Vr + 1j * Vi;           %% reconstruct V
else
    [Va, Vm] = deal(x{:});
    V = Vm .* exp(1j * Va);     %% reconstruct V
end

%% problem dimensions
nb = length(V);         %% number of buses
nl2 = length(il);       %% number of constrained lines

%% ----- evaluate constraints -----
if nl2 > 0
    flow_max = branch(il, RATE_A) / mpc.baseMVA;
    if lim_type ~= 'P'      %% typically use square of flow
        flow_max = flow_max.^2;
    end
    if lim_type == 'I'      %% current magnitude limit, |I|
        If = Yf * V;
        It = Yt * V;
        h = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
              It .* conj(It) - flow_max ];  %% branch current limits (to bus)
    else
        %% compute branch power flows
        Sf = V(branch(il, F_BUS)) .* conj(Yf * V);  %% complex power injected at "from" bus (p.u.)
        St = V(branch(il, T_BUS)) .* conj(Yt * V);  %% complex power injected at "to" bus (p.u.)
        if lim_type == '2'                      %% active power limit, P squared (Pan Wei)
            h = [ real(Sf).^2 - flow_max;       %% branch real power limits (from bus)
                  real(St).^2 - flow_max ];     %% branch real power limits (to bus)
        elseif lim_type == 'P'                  %% active power limit, P
            h = [ real(Sf) - flow_max;          %% branch real power limits (from bus)
                  real(St) - flow_max ];        %% branch real power limits (to bus
        else                                    %% apparent power limit, |S|
            h = [ Sf .* conj(Sf) - flow_max;    %% branch apparent power limits (from bus)
                  St .* conj(St) - flow_max ];  %% branch apparent power limits (to bus)
        end
    end
else
    h = zeros(0,1);
end

%%----- evaluate partials of constraints -----
if nargout > 1
    if nl2 > 0
        %% compute partials of Flows w.r.t. V
        if lim_type == 'I'                      %% current
            [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian);
        else                                    %% power
            [dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, V, mpopt.opf.v_cartesian);
        end
        if lim_type == 'P' || lim_type == '2'   %% real part of flow (active power)
            dFf_dV1 = real(dFf_dV1);
            dFf_dV2 = real(dFf_dV2);
            dFt_dV1 = real(dFt_dV1);
            dFt_dV2 = real(dFt_dV2);
            Ff = real(Ff);
            Ft = real(Ft);
        end

        if lim_type == 'P'
            %% active power
            [df_dV1, df_dV2, dt_dV1, dt_dV2] = deal(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2);
        else
            %% squared magnitude of flow (of complex power or current, or real power)
            [df_dV1, df_dV2, dt_dV1, dt_dV2] = ...
              dAbr_dV(dFf_dV1, dFf_dV2, dFt_dV1, dFt_dV2, Ff, Ft);
        end

        %% construct Jacobian of "from" and "to" branch flow ineq constraints
        dh = [ df_dV1 df_dV2;                   %% "from" flow limit
               dt_dV1 dt_dV2 ];                 %% "to" flow limit
    else
        dh = sparse(0, 2*nb);
    end
end
