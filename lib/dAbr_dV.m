function [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
						dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St)
%DABR_DV   Partial derivatives of apparent power flows w.r.t voltage.
%   [dAf_dVa, dAf_dVm, dAt_dVa, dAt_dVm] = ...
%               dAbr_dV(dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St)
%   returns four matrices containing partial derivatives of the branch
%   apparent power flows at "from" & "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If dSf_dVm is
%   a sparse matrix, the return values will be as well. The following
%   explains the expressions used to form the matrices:
%
%   Let Af refer to the apparent power at the "from" end of each line,
%   i.e. Af = abs(Sf), then ...
%
%   Partial w.r.t real power,
%       dAf/dPf = diag(real(Sf) ./ Af)
%
%   Partial w.r.t reactive power,
%       dAf/dQf = diag(imag(Sf) ./ Af)
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
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% dimensions
nl = length(Sf);

%%----- compute apparent powers -----
Af = abs(Sf);
At = abs(St);

%%----- partials w.r.t. real and reactive power flows -----
%% Careful!  Need to make partial equal to 1 for lines w/ zero flow
%%			 to avoid division by zero errors (1 comes from L'Hopital)
%% initialize to all ones
nPf = ones(nl, 1);
nQf = ones(nl, 1);
nPt = ones(nl, 1);
nQt = ones(nl, 1);
%% use actual partials for non-zero flows
i = find(Af);						%% find non-zeros of "from" flows
nPf(i) = real(Sf(i)) ./ Af(i);		
nQf(i) = imag(Sf(i)) ./ Af(i);
i = find(At);						%% find non-zeros of "to" flows
nPt(i) = real(St(i)) ./ At(i);
nQt(i) = imag(St(i)) ./ At(i);
%% put into diagonal matrices
if issparse(dSf_dVa)		%% sparse version (if dSf_dVa is sparse)
	dAf_dPf = spdiags(nPf, 0, nl, nl);
	dAf_dQf = spdiags(nQf, 0, nl, nl);
	dAt_dPt = spdiags(nPt, 0, nl, nl);
	dAt_dQt = spdiags(nQt, 0, nl, nl);
else						%% dense version
	dAf_dPf = diag(nPf);
	dAf_dQf = diag(nQf);
	dAt_dPt = diag(nPt);
	dAt_dQt = diag(nQt);
end

%% partials w.r.t. voltage magnitudes and angles
dAf_dVm = dAf_dPf * real(dSf_dVm) + dAf_dQf * imag(dSf_dVm);
dAf_dVa = dAf_dPf * real(dSf_dVa) + dAf_dQf * imag(dSf_dVa);
dAt_dVm = dAt_dPt * real(dSt_dVm) + dAt_dQt * imag(dSt_dVm);
dAt_dVa = dAt_dPt * real(dSt_dVa) + dAt_dQt * imag(dSt_dVa);

return;
