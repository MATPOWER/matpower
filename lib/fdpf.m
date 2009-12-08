function [V, converged, i] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%FDPF  Solves the power flow using a fast decoupled method.
%   [V, converged, i] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, the FDPF matrices B prime
%   and B double prime, and column vectors with the lists of bus indices
%   for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
%   vector contains the set point for generator (including ref bus)
%   buses, and the reference angle of the swing bus, as well as an initial
%   guess for remaining magnitudes and angles. mpopt is a MATPOWER options
%   vector which can be used to set the termination tolerance, maximum
%   number of iterations, and  output options (see 'help mpoption'
%   for details). Uses default options if this parameter is not given.
%   Returns the final complex voltages, a flag which indicates whether it
%   converged or not, and the number of iterations performed.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt(2);
max_it  = mpopt(4);
verbose = mpopt(31);

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
if verbose > 1
    fprintf('\niteration     max mismatch (p.u.)  ');
    fprintf('\ntype   #        P            Q     ');
    fprintf('\n---- ----  -----------  -----------');
    fprintf('\n  -  %3d   %10.3e   %10.3e', i, normP, normQ);
end
if normP < tol && normQ < tol
    converged = 1;
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% reduce B matrices
Bp = Bp([pv; pq], [pv; pq]);
Bpp = Bpp(pq, pq);

%% factor B matrices
[Lp, Up, Pp] = lu(Bp);
[Lpp, Upp, Ppp] = lu(Bpp);

%% do P and Q iterations
while (~converged && i < max_it)
    %% update iteration counter
    i = i + 1;

    %%-----  do P iteration, update Va  -----
    dVa = -( Up \  (Lp \ (Pp * P)));

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
    if verbose > 1
        fprintf('\n  P  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i-1);
        end
        break;
    end

    %%-----  do Q iteration, update Vm  -----
    dVm = -( Upp \ (Lpp \ (Ppp * Q)) );

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
    if verbose > 1
        fprintf('\n  Q  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', i, i);
        end
        break;
    end
end

if verbose
    if ~converged
        fprintf('\nFast-decoupled power flow did not converge in %d iterations.\n', i);
    end
end
