function [V, converged, i] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%FDPF  Solves the power flow using a fast decoupled method.
%   [V, CONVERGED, I] = FDPF(YBUS, SBUS, V0, BP, BPP, REF, PV, PQ, MPOPT)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, the FDPF matrices B prime
%   and B double prime, and column vectors with the lists of bus indices
%   for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
%   vector contains the set point for generator (including ref bus)
%   buses, and the reference angle of the swing bus, as well as an initial
%   guess for remaining magnitudes and angles. MPOPT is a MATPOWER options
%   vector which can be used to set the termination tolerance, maximum
%   number of iterations, and output options (see MPOPTION for details).
%   Uses default options if this parameter is not given. Returns the
%   final complex voltages, a flag which indicates whether it converged
%   or not, and the number of iterations performed.
%
%   See also RUNPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2011 by Power System Engineering Research Center (PSERC)
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

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.fd.max_it;

%% initialize
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);

%% evaluate initial mismatch
mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;
P = real(mis([pv; pq]));
Q = imag(mis(pq));

%% check tolerance
normP = norm(P, inf);
normQ = norm(Q, inf);
if mpopt.verbose > 1
    fprintf('\niteration     max mismatch (p.u.)  ');
    fprintf('\ntype   #        P            Q     ');
    fprintf('\n---- ----  -----------  -----------');
    fprintf('\n  -  %3d   %10.3e   %10.3e', i, normP, normQ);
end
if normP < tol && normQ < tol
    converged = 1;
    if mpopt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% reduce B matrices
Bp = Bp([pv; pq], [pv; pq]);
Bpp = Bpp(pq, pq);

%% factor B matrices
[Lp,  Up,  pp,  qp ] = lu(Bp,  'vector');
[Lpp, Upp, ppp, qpp] = lu(Bpp, 'vector');
[junk, iqp ] = sort(qp);
[junk, iqpp] = sort(qpp);
% [~, iqp ] = sort(qp);
% [~, iqpp] = sort(qpp);

%% do P and Q iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %%-----  do P iteration, update Va  -----
    dVa = -( Up \  (Lp \ P(pp)) );
    dVa = dVa(iqp);

    %% update voltage
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);
    if mpopt.verbose > 1
        fprintf('\n  P  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i-1);
        end
        break;
    end

    %%-----  do Q iteration, update Vm  -----
    dVm = -( Upp \ (Lpp \ Q(ppp)) );
    dVm = dVm(iqpp);

    %% update voltage
    Vm(pq) = Vm(pq) + dVm;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;
    P = real(mis([pv; pq]));
    Q = imag(mis(pq));
    
    %% check tolerance
    normP = norm(P, inf);
    normQ = norm(Q, inf);
    if mpopt.verbose > 1
        fprintf('\n  Q  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if mpopt.verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i);
        end
        break;
    end
end

if mpopt.verbose
    if ~converged
        fprintf('\nFast-decoupled power flow did not converge in %d iterations.\n', i);
    end
end
