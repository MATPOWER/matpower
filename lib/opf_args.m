function [baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
        mpopt, N, fparm, H, Cw, z0, zl, zu] = ...
    opf_args(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, ...
        mpopt, N, fparm, H, Cw, z0, zl, zu)
%OPF_ARGS  Parses and initializes OPF input arguments.
%
%   Returns the full set of initialized OPF input arguments ...
%
%   [baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt, ...
%    N, fparm, H, Cw, z0, zl, zu] = ...
%
%   ... for all of the following calling combinations ...
%
%   opf_args(mpc)
%   opf_args(mpc, mpopt)
%   opf_args(mpc, A, l, u)
%   opf_args(mpc, A, l, u, mpopt)
%   opf_args(mpc, A, l, u, mpopt, N, fparm, H, Cw)
%   opf_args(mpc, A, l, u, mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   opf_args(baseMVA, bus, gen, branch, areas, gencost)
%   opf_args(baseMVA, bus, gen, branch, areas, gencost, mpopt)
%   opf_args(baseMVA, bus, gen, branch, areas, gencost, A, l, u)
%   opf_args(baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt)
%   opf_args(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw)
%   opf_args(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   The data for the problem can be specified in one of three ways:
%   (1) a string (mpc) containing the file name of a MATPOWER case
%     which defines the data matrices baseMVA, bus, gen, branch,
%     areas and gencost.
%   (2) a struct (mpc) containing the data matrices as fields.
%   (3) the individual data matrices themselves.
%   
%   The optional user parameters for user constraints (A, l, u), user costs
%   (N, fparm, H, Cw), user variable initializer (z0), and user variable
%   limits (z0, zl, zu) can also be specified as fields in a case struct,
%   either passed in directly or defined in a case file referenced by name.
%   
%   When specified, A, l, u represent additional linear constraints on the
%   optimization variables, l <= A*[x; z] <= u. If the user specifies an A
%   matrix that has more columns than the number of "x" (OPF) variables,
%   then there are extra linearly constrained "z" variables. For an
%   explanation of the formulation used and instructions for forming the A
%   matrix, type 'help genform' or see the manual.
%
%   A generalized cost on all variables can be applied if input arguments
%   N, fparm, H and Cw are specified.  First, a linear transformation
%   of the optimization variables is defined by means of r = N * [x; z].
%   Then, to each element of r a function is applied as encoded in the
%   fparm matrix (see manual). If the resulting vector is named w,
%   then H and Cw define a quadratic cost on w: (1/2)*w'*H*w + Cw * w .
%   H and N should be sparse matrices and H should also be symmetric.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% Sort out input arguments
if isstr(baseMVA) | isstruct(baseMVA)   % passing filename or struct
  %---- opf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt, N,  fparm, H, Cw, z0, zl, zu)
  % 12  opf(casefile, Au,  lbu, ubu,    mpopt, N,       fparm, H,   Cw,  z0,    zl, zu)
  % 9   opf(casefile, Au,  lbu, ubu,    mpopt, N,       fparm, H,   Cw)
  % 5   opf(casefile, Au,  lbu, ubu,    mpopt)
  % 4   opf(casefile, Au,  lbu, ubu)
  % 2   opf(casefile, mpopt)
  % 1   opf(casefile)
  if any(nargin == [1, 2, 4, 5, 9, 12])
    casefile = baseMVA;
    if nargin == 12
      zu    = fparm;
      zl    = N;
      z0    = mpopt;
      Cw    = ubu;
      H     = lbu;
      fparm = Au;
      N     = gencost;
      mpopt = areas;
      ubu   = branch;
      lbu   = gen;
      Au    = bus;
    elseif nargin == 9
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = ubu;
      H     = lbu;
      fparm = Au;
      N     = gencost;
      mpopt = areas;
      ubu   = branch;
      lbu   = gen;
      Au    = bus;
    elseif nargin == 5
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = areas;
      ubu   = branch;
      lbu   = gen;
      Au    = bus;
    elseif nargin == 4
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = mpoption;
      ubu   = branch;
      lbu   = gen;
      Au    = bus;
    elseif nargin == 2
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = bus;
      ubu   = [];
      lbu   = [];
      Au    = sparse(0,0);
    elseif nargin == 1
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = mpoption;
      ubu   = [];
      lbu   = [];
      Au    = sparse(0,0);
    end
  else
    error('opf_args.m: Incorrect input parameter order, number or type');
  end
  mpc = loadcase(casefile);
  [baseMVA, bus, gen, branch, areas, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.areas, mpc.gencost);
  if isempty(Au) & isfield(mpc, 'A')
    [Au, lbu, ubu] = deal(mpc.A, mpc.l, mpc.u);
  end
  if isempty(N) & isfield(mpc, 'N')
    [N, fparm, H, Cw] = deal(mpc.N, mpc.fparm, mpc.H, mpc.Cw);
  end
  if isempty(z0) & isfield(mpc, 'z0')
    z0 = mpc.z0;
  end
  if isempty(zl) & isfield(mpc, 'zl')
    zl = mpc.zl;
  end
  if isempty(zu) & isfield(mpc, 'zu')
    zu = mpc.zu;
  end
else    % passing individual data matrices
  %---- opf(baseMVA, bus, gen, branch, areas, gencost, Au,   lbu, ubu, mpopt, N, fparm, H, Cw, z0, zl, zu)
  % 17  opf(baseMVA, bus, gen, branch, areas, gencost, Au,   lbu, ubu, mpopt, N, fparm, H, Cw, z0, zl, zu)
  % 14  opf(baseMVA, bus, gen, branch, areas, gencost, Au,   lbu, ubu, mpopt, N, fparm, H, Cw)
  % 10  opf(baseMVA, bus, gen, branch, areas, gencost, Au,   lbu, ubu, mpopt)
  % 9   opf(baseMVA, bus, gen, branch, areas, gencost, Au,   lbu, ubu)
  % 7   opf(baseMVA, bus, gen, branch, areas, gencost, mpopt)
  % 6   opf(baseMVA, bus, gen, branch, areas, gencost)
  if any(nargin == [6, 7, 9, 10, 14, 17])
    if nargin == 14
      zu    = [];
      zl    = [];
      z0    = [];
    elseif nargin == 10
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
    elseif nargin == 9
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = mpoption;
    elseif nargin == 7
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = Au;
      ubu   = [];
      lbu   = [];
      Au    = sparse(0,0);
    elseif nargin == 6
      zu    = [];
      zl    = [];
      z0    = [];
      Cw    = [];
      H     = [];
      fparm = [];
      N     = [];
      mpopt = mpoption;
      ubu   = [];
      lbu   = [];
      Au    = sparse(0,0);
    end
  else
    error('opf_args.m: Incorrect input parameter order, number or type');
  end
end
nw = size(N, 1);
if nw > 0
  if size(fparm, 1) ~= nw | size(H, 1) ~= nw | size(H, 2) ~= nw | ...
      length(Cw) ~= nw
    error('opf_args.m: wrong dimensions in generalized cost parameters');
  end
  if size(Au, 1) > 0 & size(N, 2) ~= size(Au, 2)
    error('opf_args.m: A and N must have the same number of columns');
  end
  % make sure N and H are sparse
  if ~issparse(N)
    N = sparse(N);
  end
  if ~issparse(H)
    H = sparse(H);
  end
end
if isempty(mpopt)
  mpopt = mpoption;
end
