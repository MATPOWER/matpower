function [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(branch, Yf, Yt, V)
%DIBR_DV   Computes partial derivatives of branch currents w.r.t. voltage.
%   [DIF_DVA, DIF_DVM, DIT_DVA, DIT_DVM, IF, IT] = DIBR_DV(BRANCH, YF, YT, V)
%   returns four matrices containing partial derivatives of the complex
%   branch currents at "from" and "to" ends of each branch w.r.t voltage
%   magnitude and voltage angle respectively (for all buses). If YF is a
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
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = ...
%           dIbr_dV(branch, Yf, Yt, V);

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define
nb = length(V);

Vnorm = V ./ abs(V);
if issparse(Yf)             %% sparse version (if Yf is sparse)
    diagV       = sparse(1:nb, 1:nb, V, nb, nb);
    diagVnorm   = sparse(1:nb, 1:nb, Vnorm, nb, nb);
else                        %% dense version
    diagV       = diag(V);
    diagVnorm   = diag(Vnorm);
end
dIf_dVa = Yf * 1j * diagV;
dIf_dVm = Yf * diagVnorm;
dIt_dVa = Yt * 1j * diagV;
dIt_dVm = Yt * diagVnorm;

%% compute currents
if nargout > 4
    If = Yf * V;
    It = Yt * V;
end
