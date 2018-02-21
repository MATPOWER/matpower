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
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
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
    [Vi, Vr] = deal(x{:});
    %% problem dimensions
    nb = length(Vi);        %% number of buses
    nl2 = length(il);       %% number of constrained lines        
    %% reconstruct V
    V = Vr + 1j*Vi;    
else
    [Va, Vm] = deal(x{:});
    %% problem dimensions
    nb = length(Va);        %% number of buses
    nl2 = length(il);       %% number of constrained lines    
    %% reconstruct V
    V = Vm .* exp(1j * Va);    
end

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
            if mpopt.opf.v_cartesian
                warning('Current magnitude limit |I| is not calculated in Cartesian coordinates')
            else
                [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V);
            end
        else                                    %% power
            if mpopt.opf.v_cartesian
                [dFf_dVi, dFf_dVr, dFt_dVi, dFt_dVr, Ff, Ft] = dSbr_dV_C(branch(il,:), Yf, Yt, V);
            else
                [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dSbr_dV_P(branch(il,:), Yf, Yt, V);
            end
        end
        if lim_type == 'P' || lim_type == '2'   %% real part of flow (active power)
            if mpopt.opf.v_cartesian
                dFf_dVi = real(dFf_dVi);
                dFf_dVr = real(dFf_dVr);
                dFt_dVi = real(dFt_dVi);
                dFt_dVr = real(dFt_dVr);
                Ff = real(Ff);
                Ft = real(Ft);                
            else
                dFf_dVa = real(dFf_dVa);
                dFf_dVm = real(dFf_dVm);
                dFt_dVa = real(dFt_dVa);
                dFt_dVm = real(dFt_dVm);
                Ff = real(Ff);
                Ft = real(Ft);                
            end
        end

        if lim_type == 'P'
            %% active power
            if mpopt.opf.v_cartesian
                [df_dVi, df_dVr, dt_dVi, dt_dVr] = deal(dFf_dVi, dFf_dVr, dFt_dVi, dFt_dVr);
            else
                [df_dVa, df_dVm, dt_dVa, dt_dVm] = deal(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm);
            end
        else
            %% squared magnitude of flow (of complex power or current, or real power)
            if mpopt.opf.v_cartesian
                [df_dVi, df_dVr, dt_dVi, dt_dVr] = ...
                  dAbr_dV_C(dFf_dVi, dFf_dVr, dFt_dVi, dFt_dVr, Ff, Ft);                
            else
                [df_dVa, df_dVm, dt_dVa, dt_dVm] = ...
                  dAbr_dV_P(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
            end
        end
        %% construct Jacobian of "from" branch flow ineq constraints        
        if mpopt.opf.v_cartesian
            dh = [ df_dVi df_dVr;                   %% "from" flow limit
                   dt_dVi dt_dVr ];                 %% "to" flow limit        
        else
            dh = [ df_dVa df_dVm;                   %% "from" flow limit
                   dt_dVa dt_dVm ];                 %% "to" flow limit            
        end
    else
        dh = sparse(0, 2*nb);
    end
end
