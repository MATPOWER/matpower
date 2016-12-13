function om = opt_model(s)
%OPT_MODEL  Constructor for optimization model class.
%   OM = OPT_MODEL
%   OM = OPT_MODEL(S)
%
%   This class implements the optimization model object used to encapsulate
%   a given optimization problem formulation. It allows for access to
%   optimization variables, constraints and costs in named blocks, keeping
%   track of the ordering and indexing of the blocks as variables,
%   constraints and costs are added to the problem.
%
%   Below are the list of available methods for use with the Opt Model class.
%   Please see the help on each individual method for more details:
%
%   Modify the OPF formulation by adding named blocks of constraints, costs
%   or variables:
%       add_constraints
%       add_costs
%       add_vars
%
%   Return the number of linear constraints, nonlinear constraints,
%   variables or cost rows, optionally for a single named block:
%       getN
%
%   Return the intial values, bounds and type for optimization variables:
%       getv
%
%   Build and return full set of linear constraints:
%       linear_constraints
%
%   Return index structure for variables, linear and nonlinear constraints
%   and costs:
%       get_idx
%
%   Build and return cost parameters and evaluate user-defined costs:
%       build_cost_params
%       get_cost_params
%       compute_cost
%
%   Save/retreive user data in the model object:
%       userdata
%
%   Display the object (called automatically when you omit the semicolon
%   at the command-line):
%       display
%
%   Return the value of an individual field:
%       get
%
%   Indentify variable, constraint or cost row indices:
%       describe_idx
%
%   The following is the structure of the data in the OPF model object.
%   Each field of .idx or .data is a struct whose field names are the names
%   of the corresponding blocks of vars, constraints or costs (found in
%   order in the corresponding .order field). The description next to these
%   fields gives the meaning of the value for each named sub-field.
%   E.g. om.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
%
%   om
%       .var        - data for optimization variable sets that make up
%                     the full optimization variable x
%           .idx
%               .i1 - starting index within x
%               .iN - ending index within x
%               .N  - number of elements in this variable set
%           .N      - total number of elements in x
%           .NS     - number of variable sets or named blocks
%           .data   - bounds and initial value data
%               .v0 - vector of initial values
%               .vl - vector of lower bounds
%               .vu - vector of upper bounds
%               .vt - scalar or vector of variable types
%                       'C' = continuous
%                       'I' = integer
%                       'B' = binary
%           .order  - struct array of names/indices for variable
%                     blocks in the order they appear in x
%               .name   - name of the block, e.g. Pg
%               .idx    - indices for name, {2,3} => Pg(2,3)
%       .nln        - data for nonlinear constraints that make up the
%                     full set of nonlinear constraints ghn(x)
%           .idx
%               .i1 - starting index within ghn(x)
%               .iN - ending index within ghn(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghn(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .order  - struct array of names/indices for nonlinear constraint
%                     blocks in the order they appear in ghn(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%       .lin        - data for linear constraints that make up the
%                     full set of linear constraints ghl(x)
%           .idx
%               .i1 - starting index within ghl(x)
%               .iN - ending index within ghl(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghl(x)
%           .NS     - number of linear constraint sets or named blocks
%           .data   - data for l <= A*xx <= u linear constraints
%               .A  - sparse linear constraint matrix
%               .l  - left hand side vector, bounding A*x below
%               .u  - right hand side vector, bounding A*x above
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for linear constraint
%                     blocks in the order they appear in ghl(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%       .cost       - data for user-defined costs
%           .idx
%               .i1 - starting row index within full N matrix
%               .iN - ending row index within full N matrix
%               .N  - number of rows in this cost block in full N matrix
%           .N      - total number of rows in full N matrix
%           .NS     - number of cost blocks
%           .data   - data for each user-defined cost block
%               .N  - see help for ADD_COSTS for details
%               .H  -               "
%               .Cw -               "
%               .dd -               "
%               .rr -               "
%               .kk -               "
%               .mm -               "
%               .vs - cell array of variable sets that define xx for this
%                     cost block, where the N for this block multiplies xx
%           .order  - struct array of names/indices for cost blocks in the
%                     order they appear in the rows of the full N matrix
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%       .userdata   - any user defined data added via USERDATA
%           .(user defined fields)

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

es = struct();
if nargin == 0
    om.var.idx.i1 = es;
    om.var.idx.iN = es;
    om.var.idx.N = es;
    om.var.N = 0;
    om.var.NS = 0;
    om.var.order = struct('name', [], 'idx', []);
    om.var.data.v0 = es;
    om.var.data.vl = es;
    om.var.data.vu = es;
    om.var.data.vt = es;

    om.nln.idx.i1 = es;
    om.nln.idx.iN = es;
    om.nln.idx.N = es;
    om.nln.N = 0;
    om.nln.NS = 0;
    om.nln.order = struct('name', [], 'idx', []);

    om.lin.idx.i1 = es;
    om.lin.idx.iN = es;
    om.lin.idx.N = es;
    om.lin.N = 0;
    om.lin.NS = 0;
    om.lin.order = struct('name', [], 'idx', []);
    om.lin.data.A = es;
    om.lin.data.l = es;
    om.lin.data.u = es;
    om.lin.data.vs = es;
    
    om.cost.idx.i1 = es;
    om.cost.idx.iN = es;
    om.cost.idx.N = es;
    om.cost.N = 0;
    om.cost.NS = 0;
    om.cost.order = struct('name', [], 'idx', []);
    om.cost.data.N = es;
    om.cost.data.H = es;
    om.cost.data.Cw = es;
    om.cost.data.dd = es;
    om.cost.data.rh = es;
    om.cost.data.kk = es;
    om.cost.data.mm = es;
    om.cost.data.vs = es;
    om.cost.params = es;
    
    om.userdata = es;

    om = class(om, 'opt_model');
elseif isa(s,'opt_model') 
    om = s;
elseif isfield(s, 'var') && isfield(s, 'nln') && ...
        isfield(s, 'lin') && isfield(s, 'cost') && isfield(s, 'userdata')
    om = class(s, 'opt_model');
else
    error('@opt_model/opt_model: input must be an OPT_MODEL object or struct with corresponding fields');
end
