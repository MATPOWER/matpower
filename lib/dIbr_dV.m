function [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V)
%DIBR_DV   Computes partial derivatives of branch currents w.r.t. voltage.
%   [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V)
%   returns four matrices containing partial derivatives of the complex
%   branch currents at "from" and "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If Yf is a
%   sparse matrix, the partial derivative matrices will be as well. Optionally
%   returns vectors containing the currents themselves. The following
%   explains the expressions used to form the matrices:
%
%   If = Yf * V;
%
%   Partials of V, Vf & If w.r.t. voltage angles
%       dV/dVa  = j * diag(V)
%       dVf/dVa = sparse(1:nl, f, j * V(f)) = j * sparse(1:nl, f, V(f))
%       dIf/dVa = Yf * dV/dVa = Yf * j * diag(V)
%
%   Partials of V, Vf & If w.r.t. voltage magnitudes
%       dV/dVm  = diag(V./abs(V))
%       dVf/dVm = sparse(1:nl, f, V(f)./abs(V(f))
%       dIf/dVm = Yf * dV/dVm = Yf * diag(V./abs(V))
%
%   Derivations for "to" bus are similar.
%

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2007 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constant
j = sqrt(-1);

%% define
nl = size(Yf, 1);
nb = length(V);

Vnorm = V ./ abs(V);
if issparse(Yf)             %% sparse version (if Yf is sparse)
    diagV       = sparse(1:nb, 1:nb, V, nb, nb);
    diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
else                        %% dense version
    diagV       = diag(V);
    diagVnorm   = diag(Vnorm);
end
dIf_dVa = Yf * j * diagV;
dIf_dVm = Yf * diagVnorm;
dIt_dVa = Yt * j * diagV;
dIt_dVm = Yt * diagVnorm;

%% compute currents
if nargout > 4
	If = Yf * V;
	It = Yt * V;
end

return;
