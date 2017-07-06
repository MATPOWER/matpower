classdef opt_model < handle
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
%   Retreive user data in the model object:
%       get_userdata
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
%       .nln        - data for (legacy) nonlinear constraints that make up the
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
%       .nle        - data for nonlinear equality constraints that make up the
%                     full set of nonlinear constraints ghne(x)
%           .idx
%               .i1 - starting index within ghne(x)
%               .iN - ending index within ghne(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghne(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .data   - data for nonlinear constraints
%               .g_dg - function handle for constraints/gradient evaluation
%               .d2G - function handle for Hessian evaluation
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for nonlinear constraint
%                     blocks in the order they appear in ghne(x)
%               .name   - name of the block, e.g. Pmis
%               .idx    - indices for name, {2,3} => Pmis(2,3)
%       .nli        - data for nonlinear inequality constraints that make up the
%                     full set of nonlinear constraints ghni(x)
%           .idx
%               .i1 - starting index within ghni(x)
%               .iN - ending index within ghni(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghni(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .data   - data for nonlinear constraints
%               .g_dg - function handle for constraints/gradient evaluation
%               .d2G - function handle for Hessian evaluation
%               .vs - cell array of variable sets that define the xx for
%                     this constraint block
%           .order  - struct array of names/indices for nonlinear constraint
%                     blocks in the order they appear in ghni(x)
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
%       .userdata   - any user defined data
%           .(user defined fields)

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%    es = struct();

    properties
        var = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'v0', struct(), ...
                'vl', struct(), ...
                'vu', struct(), ...
                'vt', struct() ) );
        nln = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ) );
        nle = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'g_dg', [], ...
                'd2G', [], ...
                'vs', struct() ) );
        nli = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'g_dg', [], ...
                'd2G', [], ...
                'vs', struct() ) );
        lin = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'A', struct(), ...
                'l', struct(), ...
                'u', struct(), ...
                'vs', struct() ) );
        cost = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'N', struct(), ...
                'H', struct(), ...
                'Cw', struct(), ...
                'dd', struct(), ...
                'rh', struct(), ...
                'kk', struct(), ...
                'mm', struct(), ...
                'vs', struct() ), ...
            'params', struct() );
        userdata = struct();
    end
    
    methods
        %% constructor
        function om = opt_model(s)
            if nargin > 0
                if isa(s, 'opt_model')
                    props = fieldnames(s);
                    for k = 1:length(props)
                        om.(props{k}) = s.(props{k});
                    end
                elseif isstruct(s)
                    props = fieldnames(om);
                    for k = 1:length(props)
                        if isfield(s, props{k})
                            om.(props{k}) = s.(props{k});
                        end
                    end
                else
                    error('@opt_model/opt_model: input must be a ''opt_model'' object or a struct');
                end
            end
        end
    end
end
