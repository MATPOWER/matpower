function P = cpf_p(parameterization, step, z, V, lam, Vprv, lamprv, pv, pq)
%CPF_P Computes partial derivatives of CPF parameterization function.
%
%   P = CPF_P(PARAMETERIZATION, STEP, Z, V, LAM, VPRV, LAMPRV, PV, PQ)
%
%   Computes the value of the parameterization function at the current
%   solution point.
%
%   Inputs:
%       PARAMETERIZATION : Value of cpf.parameterization option
%       STEP : continuation step size
%       Z : normalized tangent prediction vector from previous step
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%
%   Outputs:
%       P : value of the parameterization function at the current point
%
%   See also CPF_PREDICTOR, CPF_CORRECTOR.

%   MATPOWER
%   $Id$
%   by Shrirang Abhyankar, Argonne National Laboratory
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2013 by Power System Engineering Research Center (PSERC)
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

%% evaluate P(x0, lambda0)
if parameterization == 1        %% natural
    if lam >= lamprv
        P = lam - lamprv - step;
    else
        P = lamprv - lam - step;
    end
elseif parameterization == 2    %% arc length
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    P = sum(([Va([pv; pq]); Vm(pq); lam] - [Vaprv([pv; pq]); Vmprv(pq); lamprv]).^2) - step^2;
elseif parameterization == 3    %% pseudo arc length
    nb = length(V);
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    P = z([pv; pq; nb+pq; 2*nb+1])' * ...
        ( [Va([pv; pq]); Vm(pq); lam] - [Vaprv([pv; pq]); Vmprv(pq); lamprv] )...
        - step;
end 
