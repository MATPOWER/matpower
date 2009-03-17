function [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
                        dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St)
%DABR_DV   Partial derivatives of squared flow magnitudes w.r.t voltage.
%   [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
%               dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St)
%   returns four matrices containing partial derivatives of the square of
%   the branch flow magnitudes at "from" & "to" ends of each branch w.r.t
%   voltage magnitude and voltage angle respectively (for all buses), given
%   the flows and flow sensitivities. Flows could be complex current or
%   complex or real power. Notation is based on complex power. The following
%   explains the expressions used to form the matrices:
%
%   Let Af refer to the square of the apparent power at the "from" end of
%   each branch,
%
%       Af = abs(Sf).^2
%          = Sf .* conj(Sf)
%          = Pf.^2 + Qf.^2
%
%   then ...
%
%   Partial w.r.t real power,
%       dAf/dPf = 2 * diag(Pf)
%
%   Partial w.r.t reactive power,
%       dAf/dQf = 2 * diag(Qf)
%
%   Partial w.r.t Vm & Va
%       dAf/dVm = dAf/dPf * dPf/dVm + dAf/dQf * dQf/dVm
%       dAf/dVa = dAf/dPf * dPf/dVa + dAf/dQf * dQf/dVa
%
%   Derivations for "to" bus are similar.
%

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% dimensions
nl = length(Sf);

%%----- partials w.r.t. real and reactive power flows -----
dAf_dPf = sparse(1:nl, 1:nl, 2 * real(Sf), nl, nl);
dAf_dQf = sparse(1:nl, 1:nl, 2 * imag(Sf), nl, nl);
dAt_dPt = sparse(1:nl, 1:nl, 2 * real(St), nl, nl);
dAt_dQt = sparse(1:nl, 1:nl, 2 * imag(St), nl, nl);

%% partials w.r.t. voltage magnitudes and angles
dAf_dVm = dAf_dPf * real(dSf_dVm) + dAf_dQf * imag(dSf_dVm);
dAf_dVa = dAf_dPf * real(dSf_dVa) + dAf_dQf * imag(dSf_dVa);
dAt_dVm = dAt_dPt * real(dSt_dVm) + dAt_dQt * imag(dSt_dVm);
dAt_dVa = dAt_dPt * real(dSt_dVa) + dAt_dQt * imag(dSt_dVa);

return;
