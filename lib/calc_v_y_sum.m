function [V, Qpv, Sf, St, Sslack, iter, success] = calc_v_y_sum(Vslack,nb,nl,f,Zb,Ybf,Ybt,Yd,Sd,pv,Pg,Vg,mpopt)
% calc_v_y_sum - Solves the power flow using the admittance summation method.
% ::
%
%   [V, Qpv, Sf, St, Sslack, iter, success] = calc_v_y_sum(Vslack,nb,nl,f,Zb,Ybf,Ybt,Yd,Sd,pv,Pg,Vg,tol,iter_max)
%
%   Solves for bus voltages, generator reactive power, branch active and
%   reactive power flows and slack bus active and reactive power. The input
%   data consist of slack bus voltage, vector "from bus" indices, branch
%   impedance and shunt admittance, vector of bus shunt admittances and
%   load demand, as well as vectors with indicies of PV buses with their
%   specified voltages and active powers. It is assumed that the branches
%   are ordered using the principle of oriented ordering: indicies of
%   sending nodes are smaller then the indicies of the receiving nodes. The
%   branch index is equal to the index of their receiving node. Branch
%   addmittances are added in Yd and treated as constant admittance bus
%   loads. The applied method is admittance summation taken from:
%   Dragoslav Rajičić, Rubin Taleski, Two novel methods for radial and
%   weakly meshed network analysis, Electric Power Systems Research,
%   Volume 48, Issue 2, 15 December 1998, Pages 79-87
%   https://doi.org/10.1016/S0378-7796(98)00067-4
%
% See also radial_pf.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Mirko Todorovski
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% initialize
tol      = mpopt.pf.tol;
iter_max = mpopt.pf.radial.max_it;
vcorr    = mpopt.pf.radial.vcorr == 1;
Sd(pv) = Sd(pv) - Pg;
V = Vslack * ones(nb,1);
Vold = V;
iter = 0;
success = 0;
% ZIP load model
pw = mpopt.exp.sys_wide_zip_loads.pw;
qw = mpopt.exp.sys_wide_zip_loads.qw;
if isempty(pw)
    pw = [1 0 0];
end
if isempty(qw)
    qw = pw;
end
Sdz = real(Sd) * pw(3) + 1j * imag(Sd) * qw(3); % constant impedance
Sdi = real(Sd) * pw(2) + 1j * imag(Sd) * qw(2); % constant current
Sdp = real(Sd) * pw(1) + 1j * imag(Sd) * qw(1); % constant power
% Add artificial branch at the top of the branch list, so that the branch
% index for other branches is equal to the index of their receiving node.
 f = [0; f];
Zb = [0; Zb];
nl = nl + 1;
%% make Zpv matrix, for calculation of the PV generators reactive powers
if ~isempty(pv)
    Zpv = make_zpv(pv,nb,nl,f,Zb,Yd);
    Bpv = (imag(Zpv))^-1;
end
npv = length(pv);
Qpv = zeros(npv,1);
%% do backward-forward iterations
if mpopt.verbose > 1
    fprintf('\n it    max V mismatch (p.u.)');
    fprintf('\n----  ----------------------');
end
Ye = conj(Sdz) + Yd;
D = zeros(nl,1);
for k = nl:-1:2
    D(k) = 1 / (1 + Zb(k)*Ye(k));
    i = f(k);
    Ye(i) = Ye(i) + D(k)*Ye(k);
end
while success == 0 && iter < iter_max
    iter = iter + 1;
    % calculate load demand using actual voltages (Sdz and Yd are already
    % included in Ye)
    Vm = abs(V);
    S = Sdp + Sdi.*Vm;
    % backward sweep
    Je = conj(S./V);
    for k = nl:-1:2
        i = f(k);
        Je(i) = Je(i) + D(k)*Je(k);
    end
    % forward sweep
    for k = 2:nl
        i = f(k);
        V(k) = D(k) * (V(i)-Zb(k)*Je(k));
    end
    % check for convergence
    DU = abs(V - Vold);
    DU(isnan(DU)) = inf;
    if mpopt.verbose > 1
        fprintf('\n%3d        %10.3e', iter, max(DU));
    end
    if max(DU) > tol
        Vold = V;
        % update PV generators reactive powers
        if ~isempty(pv)
            DE = (Vg./abs(V(pv))-1).*real(V(pv)); % Rajicic (VCPF)
            DD = Bpv * DE;
            if vcorr
                DC = DD .* imag(V(pv))./real(V(pv));
                V_corr = make_vcorr(DC+1j*DD,pv,nb,nl,f,Zb);
                V = V + V_corr;
            end
            DQ = DD .* abs(V(pv)).^2 ./ real(V(pv));
            Qpv = Qpv + DQ;
            Sdp(pv) = Sdp(pv) - 1j*DQ;
        end
    else
        success = 1;
    end
end
if mpopt.verbose
    if success
        fprintf('\nAdmittance summation converged in %d iterations.\n', iter);
    else
        fprintf('\nAdmittance summation did not converge in %d iterations.\n', iter);
    end
end
%% calculate branch flows
% take out the first artificial branch
Sslack = V(1)*conj(Je(1)) + conj(Ye(1))*abs(V(1))^2;
f = f(2:end);
Zb = Zb(2:end);
% calculate branch flows
I = (V(f) - V(2:end))./Zb;
Sf = V(f).*conj(I) + conj(Ybf) .* abs(V(f)).^2;
St = V(2:end).*conj(I) - conj(Ybt) .* abs(V(2:end)).^2;