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
%   Modify the OPF formulation by adding named blocks of costs, constraints
%   or variables:
%       add_legacy_cost
%       add_quad_cost
%       add_nln_cost
%       add_lin_constraint
%       add_nln_constraint
%       add_var
%       init_indexed_name
%       add_costs (deprecated)
%       add_constraints (deprecated)
%       add_vars (deprecated)
%
%   Return the number of linear constraints, nonlinear constraints,
%   variables or cost rows, optionally for a single named block:
%       getN
%
%   Return the intial values, bounds and type for optimization variables:
%       params_var
%       getv (deprecated)
%
%   Build and return full set of linear constraints:
%       params_lin_constraint
%       linear_constraints (deprecated)
%
%   Return index structure for variables, linear and nonlinear constraints
%   and costs:
%       get_idx
%
%   Build and return cost parameters and evaluate user-defined costs:
%       params_nln_cost
%       params_quad_cost
%       params_legacy_cost
%       eval_legacy_cost
%       eval_nln_cost
%       eval_quad_cost
%       get_cost_params
%       compute_cost (deprecated)
%       build_cost_params (deprecated)
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
%       .nle        - data for nonlinear equality constraints that make up the
%                     full set of nonlinear constraints ghne(x)
%           .idx
%               .i1 - starting index within ghne(x)
%               .iN - ending index within ghne(x)
%               .N  - number of elements in this constraint set
%           .N      - total number of elements in ghne(x)
%           .NS     - number of nonlinear constraint sets or named blocks
%           .data   - data for nonlinear constraints
%               .fcn - function handle for constraints/gradient evaluation
%               .hess - function handle for Hessian evaluation
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
%               .fcn - function handle for constraints/gradient evaluation
%               .hess - function handle for Hessian evaluation
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
%       .qdc       - quadratic costs
%           .idx
%               .i1 - starting row index within quadratic costs
%               .iN - ending row index within quadratic costs
%               .N  - number of rows in this quadratic cost block
%           .N      - total number of rows in quadratic costs
%           .NS     - number of quadratic cost blocks
%           .data   - data for each quadratic cost block
%               .Q  - sparse matrix (or vector) of quadratic cost coefficients
%               .c  - column vector of linear cost coefficients
%               .k  - constant term
%               .vs - cell array of variable sets that define xx for this
%                     quadratic cost block, where sizes of Q, c, k conform to xx
%           .order  - struct array of names/indices for quadratic cost blocks
%                     in the order the were added
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%       .nlc       - general nonlinear costs
%           .idx
%               .i1 - starting row index within nonlinear costs
%               .iN - ending row index within nonlinear costs
%               .N  - number of rows in this nonlinear cost block
%                     (always equal to 1 for nonlinear cost blocks)
%           .N      - total number of rows in nonlinear costs
%           .NS     - number of nonlinear cost blocks
%           .data   - data for each nonlinear cost block
%               .fcn - function handle for cost, gradient, Hessian evaluation
%               .vs - cell array of variable sets that define xx for this
%                     nonlinear cost block, where xx is the input to the
%                     evaluation function
%           .order  - struct array of names/indices for nonlinear cost blocks
%                     in the order they were added
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%       .cost       - data for legacy user-defined costs
%           .idx
%               .i1 - starting row index within full N matrix
%               .iN - ending row index within full N matrix
%               .N  - number of rows in this cost block in full N matrix
%           .N      - total number of rows in full N matrix
%           .NS     - number of cost blocks
%           .data   - data for each user-defined cost block
%               .N  - see help for ADD_LEGACY_COST for details
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
                'fcn', [], ...
                'hess', [], ...
                'include', [], ...
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
                'fcn', [], ...
                'hess', [], ...
                'include', [], ...
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
                'vs', struct() ), ...
            'params', [] );
        qdc = struct( ...
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
                'Q', struct(), ...
                'c', struct(), ...
                'k', struct(), ...
                'vs', struct() ), ...
            'params', [] );
        nlc = struct( ...
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
                'fcn', struct(), ...
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
            'params', [] );
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
