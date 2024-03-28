function [x, f, eflag, output, lambda] = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
% qps_mips - Quadratic Program Solver based on |MIPS>|.
% ::
%
%   [x, f, exitflag, output, lambda] = ...
%       qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%
%   x = qps_mips(H, c, A, l, u)
%   x = qps_mips(H, c, A, l, u, xmin, xmax)
%   x = qps_mips(H, c, A, l, u, xmin, xmax, x0)
%   x = qps_mips(H, c, A, l, u, xmin, xmax, x0, opt)
%   x = qps_mips(problem)
%           where problem is a struct with fields:
%               H, c, A, l, u, xmin, xmax, x0, opt
%               all fields are optional, except if H is missing or all zeros,
%               then c must *not* be missing or all zeros, and linear
%               constraints or variable bounds must be supplied
%   x = qps_mips(...)
%   [x, f] = qps_mips(...)
%   [x, f, exitflag] = qps_mips(...)
%   [x, f, exitflag, output] = qps_mips(...)
%   [x, f, exitflag, output, lambda] = qps_mips(...)
%
% A wrapper function providing a standardized interface for using mips, 
% the |MIPSname| (|MIPS>|), to solve the following QP (quadratic programming)
% problem:
%
% .. math:: \min_\x f(\x) = \frac{1}{2} \trans{\x} \param{\cmat{H}} \x + \trans{\param{\rvec{c}}} \x
%   :label: eq_qps_mips_obj_fcn
%
% subject to
%
% .. math:: \param{\rvec{l}} \le \param{\cmat{A}} \x \le \param{\rvec{u}}
%   :label: eq_qps_mips_lin_constraint
% .. math:: \param{\x}_\mathrm{min} \le \x \le \param{\x}_\mathrm{max}
%   :label: eq_qps_mips_var_bounds
%
% where :eq:`eq_qps_mips_obj_fcn`--:eq:`eq_qps_mips_var_bounds` are,
% respectively, the objective function, linear constraints, and variable
% bounds.
%
% Inputs:
%   : *(all optional except that, if* ``H`` *is missing or all zeros, then*
%     ``c`` *must not be missing or all zeros, and linear constraints or
%     variable bounds must be supplied)*
%   H (double) : *(optional, default all zeros)* matrix :math:`\param{\cmat{H}}`
%       *(possibly sparse)* of quadratic cost coefficients
%   c (double) : *(optional, default all zeros)* vector of linear cost coefficients
%       :math:`\param{\rvec{c}}`
%   A, l, u (double) : *(optional, respective defaults*
%       :math:`[empty, -\infty, +\infty]` *)* matrix :math:`\param{\cmat{A}}`
%       and vectors :math:`\param{\rvec{l}}` and :math:`\param{\rvec{u}}`
%       to define the linear constraints in :eq:`eq_qps_mips_lin_constraint`
%   xmin, xmax (double) : *(optional, respective defaults*
%       :math:`[-\infty, +\infty]` *)*  lower and upper bounds,
%       :math:`\param{\x}_\mathrm{min}` and :math:`\param{\x}_\mathrm{max}`,
%       on the optimization variables :math:`\x`
%   x0 (double) : *(optional, default all zeros)* starting value of optimization vector :math:`\x`
%   opt (struct) : *(optional)* MIPS options structure (see :func:`mips` for details)
%   problem (struct) : The inputs can alternatively be supplied in a single
%       ``problem`` struct with fields corresponding to the input arguments
%       described above, namely, ``H``, ``c``, ``A``, ``l``, ``u``, ``x0``, 
%       ``xmin``, ``xmax``, and ``opt``
%
% Outputs:
%   x (double) : solution vector, :math:`\x`
%   f (double) : final objective function value, :math:`f(\x)`
%   exitflag (integer) : exit flag
%
%       - 1 = first order optimality conditions satisfied
%       - 0 = maximum number of iterations reached
%       - -1 = numerically failed
%   output (struct) : output struct with fields:
%
%       - ``iterations`` - number of iterations performed
%       - ``hist`` - struct array with trajectories of the following:
%
%         - ``feascond``, ``gradcond``, ``compcond``, ``costcond``,
%           ``gamma``, ``stepsize``, ``obj``, ``alphap``, ``alphad``
%       - ``message`` - exit message
%   lambda (struct) : struct containing the Langrange and Kuhn-Tucker
%       multipliers on the constraints, with fields:
%
%       - ``mu_l`` - lower (left-hand) limit on linear constraints
%       - ``mu_u`` - upper (right-hand) limit on linear constraints
%       - ``lower`` - lower bound on optimization variables
%       - ``upper`` - upper bound on optimization variables
%
% Note the calling syntax is almost identical to that of :func:`quadprog`
% from MathWorks' Optimization Toolbox. The main difference is that
% the linear constraints are specified with ``A``, ``l``, ``u`` instead of
% ``A``, ``B``, ``Aeq``, ``Beq``.
%
% **Example:** (problem from https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm)::
%
%   H = [   1003.1  4.3     6.3     5.9;
%           4.3     2.2     2.1     3.9;
%           6.3     2.1     3.5     4.8;
%           5.9     3.9     4.8     10  ];
%   c = zeros(4,1);
%   A = [   1       1       1       1;
%           0.17    0.11    0.10    0.18    ];
%   l = [1; 0.10];
%   u = [1; Inf];
%   xmin = zeros(4,1);
%   x0 = [1; 0; 0; 1];
%   opt = struct('verbose', 2);
%   [x, f, s, out, lambda] = qps_mips(H, c, A, l, u, xmin, [], x0, opt);
%
% See also mips.

%   MIPS
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MIPS.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mips for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(H)       %% problem struct
    p = H;
else                                %% individual args
    p = struct('H', H, 'c', c, 'A', A, 'l', l, 'u', u);
    if nargin > 5
        p.xmin = xmin;
        if nargin > 6
            p.xmax = xmax;
            if nargin > 7
                p.x0 = x0;
                if nargin > 8
                    p.opt = opt;
                end
            end
        end
    end
end

%% define nx, set default values for H and c
if ~isfield(p, 'H') || isempty(p.H) || ~any(any(p.H))
    if (~isfield(p, 'A') || isempty(p.A)) && ...
            (~isfield(p, 'xmin') || isempty(p.xmin)) && ...
            (~isfield(p, 'xmax') || isempty(p.xmax))
        error('qps_mips: LP problem must include constraints or variable bounds');
    else
        if isfield(p, 'A') && ~isempty(p.A)
            nx = size(p.A, 2);
        elseif isfield(p, 'xmin') && ~isempty(p.xmin)
            nx = length(p.xmin);
        else    % if isfield(p, 'xmax') && ~isempty(p.xmax)
            nx = length(p.xmax);
        end
    end
    p.H = sparse(nx, nx);
else
    nx = size(p.H, 1);
end
if ~isfield(p, 'c') || isempty(p.c)
    p.c = zeros(nx, 1);
end
if ~isfield(p, 'x0') || isempty(p.x0)
    p.x0 = zeros(nx, 1);
end

%%-----  run optimization  -----
p.f_fcn = @(x)qp_f(x, p.H, p.c);
[x, f, eflag, output, lambda] = mips(p);

%%-----  objective function  -----
function [f, df, d2f] = qp_f(x, H, c)
f = 0.5 * x' * H * x + c' * x;
if nargout > 1
    df = H * x + c;
    if nargout > 2
        d2f = H;
    end
end
