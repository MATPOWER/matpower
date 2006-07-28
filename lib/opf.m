function [buso, gen, branch, f, success, info, et, g, jac] = opf(baseMVA, bus,...
          gen, branch, areas, gencost, Au, lbu, ubu, mpopt, N, fparm, H, Cw)
%OPF  Solves an optimal power flow.
%
%   For an AC OPF, if the OPF algorithm is not set explicitly in the options,
%   it will choose the best available solver, searching in the following order:
%   MINOPF, fmincon, constr.
%
%   [bus, gen, branch, f, success] = opf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = opf(casefile, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, areas, ...
%                                    gencost, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, areas, ...
%                                    gencost, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, areas, ...
%                                    gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw)
%
%   [bus, gen, branch, f, success, info, et, g, jac] = opf(casefile)
%
%   The data for the problem can be specified in one of 3 ways: (1) the name of
%   a case file which defines the data matrices baseMVA, bus, gen, branch,
%   areas and gencost, (2) a struct containing the data matrices as fields, or
%   (3) the data matrices themselves.
%
%   When specified, A, l, u represent additional linear constraints on the
%   optimization variables, l <= A*[x; z] <= u. For an explanation of the
%   formulation used and instructions for forming the A matrix, type
%   'help genform'.
%
%   A generalized cost on all variables can be applied if input arguments
%   N, fparm, H and Cw are specified.  First, a linear transformation
%   of the optimization variables is defined by means of r = N * [x; z].
%   Then, to each element of r a function is applied as encoded in the
%   fparm matrix (see manual or type 'help generalcost').  If the
%   resulting vector is now named w, then H and Cw define a quadratic
%   cost on w:  (1/2)*w'*H*w + Cw * w . H and N should be sparse matrices
%   and H should also be symmetric.
%
%   The additional linear constraints and generalized cost are only available
%   for solvers which use the generalized formulation, namely fmincon and
%   MINOPF.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the constraint vector (g) and the
%   Jacobian matrix (jac).
%
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the number of "x" (OPF) variables, then there are extra linearly
%   constrained "z" variables.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2006 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% Sort out input arguments
if isstr(baseMVA) | isstruct(baseMVA)   % passing filename or struct
  % 14  opf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt, N, fparm, H, Cw)
  % 9   opf(casefile, Au,  lbu, ubu,    mpopt, N,       fparm, H,   Cw)
  % 5   opf(casefile, Au,  lbu, ubu,    mpopt)
  % 4   opf(casefile, Au,  lbu, ubu)
  % 2   opf(casefile, mpopt)
  % 1   opf(casefile)
  if any(nargin == [1, 2, 4, 5, 9])
    casefile = baseMVA;
    if nargin == 9
      N     = gencost;
      fparm = Au;
      H     = lbu;
      Cw    = ubu;
    else
      N     = [];
      fparm = [];
      H     = [];
      Cw    = [];
    end
    if nargin < 4
      Au  = sparse(0,0);
      lbu = [];
      ubu = [];
    else
      Au  = bus;
      lbu = gen;
      ubu = branch;
    end
    if nargin == 9 | nargin == 5
      mpopt = areas;
    elseif nargin == 2
      mpopt = bus;
    else
      mpopt = [];
    end
  else
    error('opf.m: Incorrect input parameter order, number or type');
  end
  [baseMVA, bus, gen, branch, areas, gencost] = loadcase(casefile);
else    % passing individual data matrices
  % 14  opf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt, N, fparm, H, Cw)
  % 10  opf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt)
  % 9   opf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu)
  % 7   opf(baseMVA,  bus, gen, branch, areas, gencost, mpopt)
  % 6   opf(baseMVA,  bus, gen, branch, areas, gencost)
  if any(nargin == [6, 7, 9, 10, 14])
    if nargin ~= 14
      N     = [];
      fparm = [];
      H     = [];
      Cw    = [];
    end
    if nargin == 7
      mpopt = Au;
    elseif nargin == 6 | nargin == 9
      mpopt = [];
    end
    if nargin < 9
      Au  = sparse(0,0);
      lbu = [];
      ubu = [];
    end
  else
    error('opf.m: Incorrect input parameter order, number or type');
  end
end
if size(N, 1) > 0
  if size(N, 1) ~= size(fparm, 1) | size(N, 1) ~= size(H, 1) | ...
     size(N, 1) ~= size(H, 2) | size(N, 1) ~= length(Cw)
    error('opf.m: wrong dimensions in generalized cost parameters');
  end
  if size(Au, 1) > 0 & size(N, 2) ~= size(Au, 2)
    error('opf.m: A and N must have the same number of columns');
  end
end
if isempty(mpopt)
  mpopt = mpoption;
end

%%----- initialization -----
%% options
verbose = mpopt(31);
npts = mpopt(14);       %% number of points to evaluate when converting
                        %% polynomials to piece-wise linear

%% define constants
j = sqrt(-1);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%%-----  check/convert costs, set default algorithm  -----
%% get cost model, check consistency
model = gencost(:, MODEL);
comgen = find(gen(:, GEN_STATUS) > 0);
if size(gencost, 1) == 2*size(gen,1)
  comgen = [comgen; comgen];
end
i_pwln = find(model(comgen) == PW_LINEAR);
i_poly = find(model(comgen) == POLYNOMIAL);

% Start clock
t1 = clock;

%% set algorithm
dc = mpopt(10);
if dc % DC OPF
  [bus, gen, branch, f, success, info, et] = dcopf(baseMVA, bus, gen, ...
                                            branch, areas, gencost, mpopt);
  g = []; jac = [];     %% not currently available from DC OPF
else % AC optimal power flow requested
  if any(model ~= PW_LINEAR & model ~= POLYNOMIAL)
    error('opf.m: unknown generator cost model in gencost data');
  end
  if mpopt(11) == 0  % OPF_ALG not set, choose best option
    if have_fcn('minopf')
      mpopt(11) = 500; % MINOS generalized
    elseif have_fcn('fmincon')
      mpopt(11) = 520; % FMINCON generalized
    %% use default for this cost model
    elseif any(i_pwln)      %% some piece-wise linear, use appropriate alg
      mpopt(11) = mpopt(13);
      if any(i_poly) & verbose
        fprintf('opf.m: not all generators use same cost model, all will be converted\n       to piece-wise linear\n');
      end
    else                    %% must all be polynomial
        mpopt(11) = mpopt(12);
    end
  end
  alg = mpopt(11);
  formulation = opf_form(alg); % 1, 2 or 5

  %% check cost model/algorithm consistency
  if any( i_pwln ) & formulation == 1
    error(sprintf('opf.m: algorithm %d does not handle piece-wise linear cost functions', alg));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Eventually put code here to fit polynomials to piece-wise linear as needed.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end

  if (formulation ~= 5) & ~isempty(Au)
    error('opf.m: Selected algorithm cannot handle general linear constraints');
  end
  
  %% convert polynomials to piece-wise linear
  if any(i_poly)  & formulation == 2
    if verbose
        fprintf('converting from polynomial to piece-wise linear cost model\n');
    end
    [pcost, qcost] = pqcost(gencost, size(gen, 1));
    i_poly = find(pcost(:, MODEL) == POLYNOMIAL);
    tmp = poly2pwl(pcost(i_poly, :), gen(i_poly, PMIN), gen(i_poly, PMAX), npts);
    pcost(i_poly, 1:size(tmp,2)) = tmp;
    if ~isempty(qcost)
        i_poly = find(qcost(:, MODEL) == POLYNOMIAL);
        tmp = poly2pwl(qcost(i_poly, :), gen(i_poly, QMIN), gen(i_poly, QMAX), npts);
        qcost(i_poly, 1:size(tmp,2)) = tmp;
    end
    gencost = [pcost; qcost];
  end

  %%-----  run opf  -----
  if formulation == 5 % Generalized
    if alg == 500       % MINOS
      if ~have_fcn('minopf')
        error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
            'MINOPF (see http://www.pserc.cornell.edu/minopf/)']);
      end
      [bus, gen, branch, f, success, info, et, g, jac] = mopf(baseMVA, ...
          bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, N, fparm, H, Cw);
    elseif alg == 520   % FMINCON
      if ~have_fcn('fmincon')
        error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
            'fmincon (Optimization Toolbox 2.x or later)']);
      end
      [bus, gen, branch, f, success, info, et, g, jac] = fmincopf(baseMVA, ...
          bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, N, fparm, H, Cw);
    end
  else
    if opf_slvr(alg) == 0           %% use CONSTR
      if ~have_fcn('constr')
        error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
            'constr (Optimization Toolbox 1.x/2.x)']);
      end
      %% set some options
      if mpopt(19) == 0
        mpopt(19) = 2 * size(bus,1) + 150;  %% set max number of iterations for constr
      end
   
      %% run optimization
      [bus, gen, branch, f, success, info, et, g, jac] = copf(baseMVA, ...
              bus, gen, branch, areas, gencost, mpopt);
  
    else                            %% use LPCONSTR
      [bus, gen, branch, f, success, info, et, g, jac] = lpopf(baseMVA, ...
              bus, gen ,branch, areas, gencost, mpopt);
    end
  end
end
    
%% compute elapsed time
et = etime(clock, t1);
if (nargout == 0) & ( success )
  printpf(baseMVA, bus, gen, branch, f, success, et, 1, mpopt);
end

if nargout
  buso = bus;
end

return;
