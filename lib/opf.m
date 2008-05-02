function [buso, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
    opf(varargin)
%OPF  Solves an optimal power flow.
%
%   For an AC OPF, if the OPF algorithm is not set explicitly in the options,
%   it will choose the best available solver, searching in the following order:
%   minopf, pdipmopf, fmincon, constr.
%
%   [bus, gen, branch, f, success] = opf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = opf(casefile, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw)
%
%   [bus, gen, branch, f, success] = opf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw, z0, zl, zu)
%
%   [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = opf(casefile)
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
%   fparm matrix (see manual).  If the resulting vector is now named w,
%   then H and Cw define a quadratic cost on w: (1/2)*w'*H*w + Cw * w .
%   H and N should be sparse matrices and H should also be symmetric.
%
%   The additional linear constraints and generalized cost are only available
%   for solvers which use the generalized formulation, namely fmincon,
%   MINOPF and TSPOPF.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the constraint vector (g), the
%   Jacobian matrix (jac), and the vector of variables (xr) as well 
%   as the constraint multipliers (pimul).
%
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the number of "x" (OPF) variables, then there are extra linearly
%   constrained "z" variables.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% process input arguments
[baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu] = opf_args(varargin{:});

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

%% initialize optional output args
g = []; jac = []; xr = []; pimul = [];

% Start clock
t1 = clock;

%% set algorithm
dc = mpopt(10);
if dc % DC OPF
  nb = size(bus, 1);
  ng = size(gen, 1);
  if size(Au, 1) > 0
    %% reduce columns of Au and N, assuming they were for the AC problem
    %% the AC problem unless there aren't enough columns
    nxz = size(Au, 2);
    if nxz >= 2*nb + 2*ng
      %% define indexing of optimization variable vector
      k = 0;
      thbas   = k + 1;    k = k + nb;     thend = k;      %% voltage angles
      vbas    = k + 1;    k = k + nb;     vend  = k;      %% voltage magnitudes
      pgbas   = k + 1;    k = k + ng;     pgend = k;      %% real power injections
      qgbas   = k + 1;    k = k + ng;     qgend = k;      %% reactive power injections
      %% make sure there aren't any constraints or costs on V or Qg
      if any(any(Au(:, [vbas:vend qgbas:qgend])))
        error('opf: Attempting to solve DC OPF with user constraints on V or Qg');
      end
      dcc = [thbas:thend pgbas:pgend];
      if nxz > qgend
        dcc = [dcc (qgend+1):nxz];
      end
      Au = Au(:, dcc);
      nw = size(N, 1);
      if nw > 0
        if any(any(N(:, [vbas:vend qgbas:qgend])))
          error('opf: Attempting to solve DC OPF with user costs on V or Qg');
        end
        N = N(:, dcc);
      end
    end
  end
  if nargout > 7
    [bus, gen, branch, f, success, info, et, xr, pimul] = ...
        dcopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
            mpopt, N, fparm, H, Cw, z0, zl, zu);
  else
    [bus, gen, branch, f, success, info, et] = ...
        dcopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
            mpopt, N, fparm, H, Cw, z0, zl, zu);
  end
else % AC optimal power flow requested
  if any(model ~= PW_LINEAR & model ~= POLYNOMIAL)
    error('opf.m: unknown generator cost model in gencost data');
  end
  if mpopt(11) == 0  % OPF_ALG not set, choose best option
    if have_fcn('minopf')
      mpopt(11) = 500; % MINOS generalized
    elseif have_fcn('pdipmopf')
      mpopt(11) = 540; % PDIPM generalized
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
      if nargout > 7
        [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
            mopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
                mpopt, N, fparm, H, Cw, z0, zl, zu);
      else
        [bus, gen, branch, f, success, info, et] = mopf(baseMVA, ...
            bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
            N, fparm, H, Cw, z0, zl, zu);
      end
    elseif alg == 520   % FMINCON
      if ~have_fcn('fmincon')
        error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
            'fmincon (Optimization Toolbox 2.x or later)']);
      end
      mlver = ver('matlab');
      if str2num(mlver.Version(1)) < 7    %% anonymous functions not available
        fmc = @fmincopf6;
      else
        fmc = @fmincopf;
      end
      if nargout > 7
        [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
            feval(fmc, baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
                mpopt, N, fparm, H, Cw, z0, zl, zu);
      else
        [bus, gen, branch, f, success, info, et] = feval(fmc, baseMVA, ...
            bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
            N, fparm, H, Cw, z0, zl, zu);
      end
    elseif alg == 540 | alg == 545 | alg == 550  % PDIPM_OPF, SCPDIPM_OPF, or TRALM_OPF
      if alg == 540       % PDIPM_OPF
        if ~have_fcn('pdipmopf')
          error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
              'PDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)']);
        end
      elseif alg == 545       % SCPDIPM_OPF
        if ~have_fcn('scpdipmopf')
          error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
              'SCPDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)']);
        end
      elseif alg == 550       % TRALM_OPF
        if ~have_fcn('tralmopf')
          error(['opf.m: OPF_ALG ', num2str(alg), ' requires ', ...
              'TRALMOPF (see http://www.pserc.cornell.edu/tspopf/)']);
        end
      end
      if nargout > 7
        [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = ...
            tspopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
                mpopt, N, fparm, H, Cw, z0, zl, zu);
      else
        [bus, gen, branch, f, success, info, et] = tspopf(baseMVA, ...
            bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
            N, fparm, H, Cw, z0, zl, zu);
      end
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
      if nargout > 8
        [bus, gen, branch, f, success, info, et, g, jac] = ...
            copf(baseMVA, bus, gen, branch, areas, gencost, mpopt);
      else
        [bus, gen, branch, f, success, info, et, g] = ...
            copf(baseMVA, bus, gen, branch, areas, gencost, mpopt);
      end
    else                            %% use LPCONSTR
      [bus, gen, branch, f, success, info, et] = lpopf(baseMVA, ...
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
