function [dP_dV, dP_dlam] = cpf_p_jac(parameterization, z, V, lam, Vprv, lamprv, pv, pq)
%CPF_P_JAC Computes partial derivatives of CPF parameterization function.
%
%   [DP_DV, DP_DLAM ] = CPF_P_JAC(PARAMETERIZATION, Z, V, LAM, ...
%                                                   VPRV, LAMPRV, PV, PQ)
%
%   Computes the partial derivatives of the continuation power flow
%   parameterization function w.r.t. bus voltages and the continuation
%   parameter lambda.
%
%   Inputs:
%       PARAMETERIZATION : Value of cpf.parameterization option.
%       Z : normalized tangent prediction vector from previous step
%       V : complex bus voltage vector at current solution
%       LAM : scalar lambda value at current solution
%       VPRV : complex bus voltage vector at previous solution
%       LAMPRV : scalar lambda value at previous solution
%       PV : vector of indices of PV buses
%       PQ : vector of indices of PQ buses
%
%   Outputs:
%       DP_DV : partial of parameterization function w.r.t. voltages
%       DP_DLAM : partial of parameterization function w.r.t. lambda
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

if parameterization == 1        %% natural
    npv = length(pv);
    npq = length(pq);
    dP_dV = zeros(1, npv+2*npq);
    if lam >= lamprv
        dP_dlam = 1.0;
    else
        dP_dlam = -1.0;
    end
elseif parameterization == 2    %% arc length
    Va = angle(V);
    Vm = abs(V);
    Vaprv = angle(Vprv);
    Vmprv = abs(Vprv);
    dP_dV = 2*([Va([pv; pq]); Vm(pq)] - [Vaprv([pv; pq]); Vmprv(pq)])';
    if lam == lamprv    %% first step
        dP_dlam = 1.0;  %% avoid singular Jacobian that would result
                        %% from [dP_dV, dP_dlam] = 0
    else
        dP_dlam = 2*(lam-lamprv);
    end
elseif parameterization == 3    %% pseudo arc length
    nb = length(V);
    dP_dV = z([pv; pq; nb+pq])';
    dP_dlam = z(2*nb+1);
end 
