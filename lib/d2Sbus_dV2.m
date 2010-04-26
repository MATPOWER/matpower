function [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam)
%D2SBUS_DV2   Computes 2nd derivatives of power injection w.r.t. voltage.
%   [GAA, GAV, GVA, GVV] = D2SBUS_DV2(YBUS, V, LAM) returns 4 matrices
%   containing the partial derivatives w.r.t. voltage angle and magnitude
%   of the product of a vector LAM with the 1st partial derivatives of the
%   complex bus power injections. Takes sparse bus admittance matrix YBUS,
%   voltage vector V and nb x 1 vector of multipliers LAM. Output matrices
%   are sparse.
%
%   Example:
%       [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
%       [Gaa, Gav, Gva, Gvv] = d2Sbus_dV2(Ybus, V, lam);
%
%   Here the output matrices correspond to:
%       Gaa = (d/dVa (dSbus_dVa.')) * lam
%       Gav = (d/dVm (dSbus_dVa.')) * lam
%       Gva = (d/dVa (dSbus_dVm.')) * lam
%       Gvv = (d/dVm (dSbus_dVm.')) * lam

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
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

n = length(V);
Ibus    = Ybus * V;
diaglam = sparse(1:n, 1:n, lam, n, n);
diagV   = sparse(1:n, 1:n, V, n, n);

A = sparse(1:n, 1:n, lam .* V, n, n);
B = Ybus * diagV;
C = A * conj(B);
D = Ybus' * diagV;
E = conj(diagV) * (D * diaglam - sparse(1:n, 1:n, D*lam, n, n));
F = C - A * sparse(1:n, 1:n, conj(Ibus), n, n);
G = sparse(1:n, 1:n, ones(n, 1)./abs(V), n, n);

Gaa = E + F;
Gva = 1j * G * (E - F);
Gav = Gva.';
Gvv = G * (C + C.') * G;
