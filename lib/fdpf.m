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
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt.pf.tol;
max_it  = mpopt.pf.fd.max_it;
lu_vec  = have_fcn('lu_vec');

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
mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
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
if lu_vec
    [Lp,  Up,  pp,  qp ] = lu(Bp,  'vector');
    [Lpp, Upp, ppp, qpp] = lu(Bpp, 'vector');
    [junk, iqp ] = sort(qp);
    [junk, iqpp] = sort(qpp);
    % [~, iqp ] = sort(qp);
    % [~, iqpp] = sort(qpp);
else
    [Lp, Up, Pp] = lu(Bp);
    [Lpp, Upp, Ppp] = lu(Bpp);
end

%% do P and Q iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %%-----  do P iteration, update Va  -----
    if lu_vec
        dVa = -( Up \  (Lp \ P(pp)) );
        dVa = dVa(iqp);
    else
        dVa = -( Up \  (Lp \ (Pp * P)));
    end

    %% update voltage
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
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
    if lu_vec
        dVm = -( Upp \ (Lpp \ Q(ppp)) );
        dVm = dVm(iqpp);
    else
        dVm = -( Upp \ (Lpp \ (Ppp * Q)) );
    end

    %% update voltage
    Vm(pq) = Vm(pq) + dVm;
    V = Vm .* exp(1j * Va);

    %% evalute mismatch
    mis = (V .* conj(Ybus * V) - Sbus(Vm)) ./ Vm;
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
