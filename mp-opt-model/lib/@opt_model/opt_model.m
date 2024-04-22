classdef opt_model < mp_idx_manager
% opt_model - MP-Opt-Model optimization model class.
% ::
%
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
%   Modify the OPF formulation by adding named blocks of costs, constraints,
%   or variables:
%       add_quad_cost
%       add_nln_cost
%       add_lin_constraint
%       add_nln_constraint
%       add_var
%       init_indexed_name
%
%   Modify parameters of existing costs, constraints, or variables:
%       set_params
%
%   Return the number of linear constraints, nonlinear constraints, or
%   variables, optionally for a single named block:
%       getN
%
%   Return index structure for variables, linear and nonlinear constraints
%   and costs:
%       get_idx
%
%   Return the intial values, bounds and type for optimization variables:
%       params_var
%
%   Build and return constraint parameters and evaluate constraints:
%       params_lin_constraint
%       params_nln_constraint
%       eval_nln_constraint
%       eval_nln_constraint_hess
%
%   Build and return cost parameters and evaluate costs:
%       params_nln_cost
%       params_quad_cost
%       eval_nln_cost
%       eval_quad_cost
%
%   Determine model type:
%       problem_type
%
%   Solve the model, access, and display the solution:
%       solve
%       display_soln
%       get_soln
%       parse_soln
%       is_solved
%       has_parsed_soln
%
%   Retreive user data in the model object:
%       get_userdata
%
%   Display the object (called automatically when you omit the semicolon
%   at the command-line):
%       display
%
%   Utility methods for dealing with varsets:
%       varsets_cell2struct
%       varsets_idx
%       varsets_len
%       varsets_x
%
%   Return the value of an individual field:
%       get
%
%   Indentify variable, constraint or cost row indices:
%       describe_idx
%
%   Make a shallow copy of the object:
%       copy
%
%   The following is the structure of the data in the Opt-Model object.
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
%           .params - cache for previously assembled aggregate parameters
%               .v0  - aggregate vector of variable initial values
%               .vl  - aggregate vector of variable lower bounds
%               .vu  - aggregate vector of variable upper bounds
%               .vt  - aggregate vector of variable types
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
%           .params - cache for previously assembled aggregate parameters
%               .A  - aggregate sparse linear constraint matrix
%               .l  - aggregate left hand side vector, bounding A*x below
%               .u  - aggregate right hand side vector, bounding A*x above
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
%           .params - cache for previously assembled aggregate parameters
%               .Q  - aggregate sparse matrix of quadratic cost coefficients
%               .c  - aggregate column vector of linear cost coefficients
%               .k  - aggregate constant term
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
%       .prob_type  - used to cache the return value of problem_type()
%       .soln       - struct containing the output of the solve() method
%                     with the following fields
%           .eflag  - exit flag
%           .output - output struct with the following fields
%               .alg     - algorithm code of solver used
%               (others) - solver specific fields
%           .x      - solution vector
%           .f      - final (objective) function value
%           .jac    - final Jacobian matrix (if available, for LEQ/NLEQ probs)
%           .lambda - Lagrange and Kuhn-Tucker multipliers on constraints
%       .userdata   - any user defined data
%           .(user defined fields)
%
% See also mp_idx_manager.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        var             %% variables
        lin             %% linear constraints
        nle             %% nonlinear equality constraints
        nli             %% nonlinear inequality constraints
        qdc             %% quadratic costs
        nlc             %% general nonlinear costs
        prob_type = ''; %% problem type
        
        %% results of solve()
        soln = struct( ...
            'eflag', [], ...    %% exit flag
            'output', [], ...   %% algorithm code & solver-specific fields
            'x', [], ...        %% solution vector
            'f', [], ...        %% final (objective) function value
            'jac', [], ...      %% Jacobian (if available) for LEQ/NLEQ
            'lambda', [] );     %% constraint shadow prices
    end     %% properties

    methods
        function om = opt_model(varargin)
            % Constructor.
            % ::
            %
            %   om = opt_model()
            %   om = opt_model(a_struct)
            %   om = opt_model(an_om)

            %% call parent constructor
            om@mp_idx_manager(varargin{:});

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if om.var is empty) in ADD_VAR(), DISPLAY() and
            %%              INIT_INDEXED_NAME(), after object construction,
            %%              but before object use.
        end

        function om = def_set_types(om)
            % Define set types ``var``, ``lin``, ``nle``, ``nli``, ``qdc``, ``nlc``.

            om.set_types = struct(...
                    'var', 'VARIABLES', ...
                    'lin', 'LINEAR CONSTRAINTS', ...
                    'nle', 'NONLIN EQ CONSTRAINTS', ...
                    'nli', 'NONLIN INEQ CONSTRAINTS', ...
                    'qdc', 'QUADRATIC COSTS', ...
                    'nlc', 'GEN NONLIN COSTS' ...
                );
        end

        function om = init_set_types(om)
            % Initialize data structures for each set type.

            %% call parent to create base data structures for each type
            init_set_types@mp_idx_manager(om);

            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            om.var.data = struct( ...
                'v0', es, ...
                'vl', es, ...
                'vu', es, ...
                'vt', es );
            om.var.params = [];
            om.nle.data = struct( ...
                'fcn', [], ...
                'hess', [], ...
                'include', [], ...
                'vs', es );
            om.nli.data = struct( ...
                'fcn', [], ...
                'hess', [], ...
                'include', [], ...
                'vs', es );
            om.lin.data = struct( ...
                'A', es, ...
                'l', es, ...
                'u', es, ...
                'tr', es, ...
                'vs', es );
            om.lin.params = [];
            om.qdc.data = struct( ...
                'Q', es, ...
                'c', es, ...
                'k', es, ...
                'vs', es );
            om.qdc.params = [];
            om.nlc.data = struct( ...
                'fcn', es, ...
                'vs', es );
        end

        function TorF = is_solved(om)
            % Return true if model has been solved.

            TorF = ~isempty(om.soln.eflag);
        end

        function TorF = has_parsed_soln(om)
            % Return true if model has a parsed solution.

            TorF = isfield(om.soln, 'var');
        end
    end     %% methods
end         %% classdef
