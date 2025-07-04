classdef opt_model < handle
% mp.opt_model - Mathematical programming and optimization model class.
% ::
%
%   mm = mp.opt_model
%   mm = mp.opt_model(s)
%   mm = mp.opt_model(mm0)
%   mm = mp.opt_model(om)
%
% This class implements the model object used to encapsulate a given
% mathematical programming or optimization problem formulation. It allows
% for access to variables, constraints and costs in named blocks,
% keeping track of the ordering and indexing of the blocks as variables,
% constraints and costs are added to the model.
%
% mp.opt_model Properties:
%   * var - MP Set Manager for variables
%   * lin - MP Set Manager for linear constraints
%   * qcn - MP Set Manager for quadratic constraints
%   * nle - MP Set Manager for nonlinear equality constraints
%   * nli - MP Set Manager for nonlinear inequality constraints
%   * qdc - MP Set Manager for quadratic costs
%   * nlc - MP Set Manager for general nonlinear costs
%   * prob_type - problem type, cached result of problem_type()
%   * soln - results of solve()
%   * userdata - arbitrary user data
%
% mp.opt_model Methods:
%   * opt_model - constructor, optionally copying from struct or existing object
%   * get_set_types - list of names of properties of set types managed by this class
%   * copy - make a duplicate of the object
%   * to_struct - convert object data *to* a struct
%   * from_struct - copy object data *from* a struct
%   * get_idx - return ``idx`` struct for vars, constraints, costs
%   * get_user_data - used to retrieve values of user data
%   * problem_type - return string identifying type of mathematical program
%   * is_mixed_integer - return true if model is mixed integer, false otherwise
%   * is_solved - return true if model has been solved
%   * solve - solve the model
%   * has_parsed_soln - return true if model has a parsed solution
%   * parse_soln - parse solution vector and shadow prices by all named sets
%   * display - display the object
%   * display_header - display header, before each set type
%   * display_footer - display footer, after each set type
%   * display_soln - display solution values
%
% Example::
%
%   mm = mp.opt_model;
%   mm.var.add('y', 2, y0, ymin);
%   mm.var.add('z', 2, z0, [], zmax);
%   mm.lin.add(mm.var, 'lincon1', A1, b1, b1);
%   mm.lin.add(mm.var, 'lincon2', A2, [], u2, {'y'});
%   mm.qdc.add(mm.var, 'cost', H, []);
%
%   prob_type = mm.problem_type();
%
%   if ~mm.is_solved()
%       [x, f, exitflag, output, lambda] = mm.solve();
%   end
%   mm.display_soln();
%
% Replaces the now deprecated legacy opt_model class and its base class,
% mp_idx_manager.
%
% See also mp.set_manager.

%   MP-Opt-Model
%   Copyright (c) 2008-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

    properties
        var             % (mp.sm_variable) MP Set Manager for variables
        lin             % (mp.sm_lin_constraint) MP Set Manager for linear constraints
        qcn             % (mp.sm_quad_constraint) MP Set Manager for quadratic constraints
        nle             % (mp.sm_nln_constraint) MP Set Manager for nonlinear equality constraints
        nli             % (mp.sm_nln_constraint) MP Set Manager for nonlinear inequality constraints
        qdc             % (mp.sm_quad_cost) MP Set Manager for quadratic costs
        nlc             % (mp.sm_nln_cost) MP Set Manager for general nonlinear costs
        prob_type = ''; % *(char array)* problem type, cached result of problem_type()

        % *(struct)* results of solve(), with fields:
        %
        %   - ``'eflag'`` - exit flag
        %   - ``'output'`` - algorithm code (``'alg'``) & solver-specific fields
        %   - ``'x'`` - solution vector
        %   - ``'f'`` - final (objective) function value
        %   - ``'jac'`` - Jacobian (if available) for LEQ/NLEQ
        %   - ``'lambda'`` - struct of constraint shadow prices
        soln = struct( ...
            'eflag', [], ...
            'output', [], ...
            'x', [], ...
            'f', [], ...
            'jac', [], ...
            'lambda', [] );
        userdata = struct();    % *(struct)* arbitrary user data
    end     %% properties

    methods
        function mm = opt_model(varargin)
            % Constructor.
            % ::
            %
            %   mm = mp.opt_model()
            %   mm = mp.opt_model(a_struct)
            %   mm = mp.opt_model(an_mm)
            %   mm = mp.opt_model(an_om)
            %
            % Optionally copies data from (duplicates) a struct or existing
            % mathematical model object.
            %
            % Input:
            %   a_struct (struct) : a struct with the stucture obtained by
            %       converting an mp.opt_model object and its sub-objects to a
            %       struct
            %   an_mm (mp.opt_model) : existing mathematical model
            %       (mp.opt_model) object
            %   an_om (:class:`opt_model`) : existing *legacy* mathematical
            %       model (:class:`opt_model`) object

            if nargin > 0
                s = varargin{1};
                if isa(s, 'mp.opt_model') || isa(s, 'opt_model')
                    if have_feature('octave')
                        s1 = warning('query', 'Octave:classdef-to-struct');
                        warning('off', 'Octave:classdef-to-struct');
                    end
                    props = fieldnames(s);
                    if have_feature('octave')
                        warning(s1.state, 'Octave:classdef-to-struct');
                    end
                    [~, k] = ismember('set_types', props);
                    if k
                        props(k) = [];  %% remove 'set_types'
                    end
                    for k = 1:length(props)
                        mm = copy_prop(s, mm, props{k});
                    end
                elseif isstruct(s)
                    props = fieldnames(mm);
                    for k = 1:length(props)
                        if isfield(s, props{k})
                            mm = copy_prop(s, mm, props{k});
                        end
                    end
                    st = mm.get_set_types();
                    for k = 1:length(st)
                        if isstruct(mm.(st{k}))
                            mm.(st{k}) = s.(st{k}).to_struct();
                        end
                    end
                else
                    error('mp.opt_model: input must be an ''mp.opt_model'' object or a struct');
                end
            end

            if isempty(mm.var)      %% skip for copy constructor
                mm.var = mp.sm_variable('VARIABLES');
                mm.lin = mp.sm_lin_constraint('LINEAR CONSTRAINTS');
                mm.qcn = mp.sm_quad_constraint('QUADRATIC CONSTRAINTS');
                mm.nle = mp.sm_nln_constraint('NONLIN EQ CONSTRAINTS');
                mm.nli = mp.sm_nln_constraint('NONLIN INEQ CONSTRAINTS');
                mm.qdc = mp.sm_quad_cost('QUADRATIC COSTS');
                mm.nlc = mp.sm_nln_cost('GEN NONLIN COSTS');
            end
        end

        function st = get_set_types(mm)
            % List of names of properties of set types managed by this class.
            % ::
            %
            %   st = mm.get_set_types();
            %
            % Output:
            %   st (cell array) : list of set types, namely:
            %       ``{'var', 'lin', 'qcn', 'nle', 'nli', 'qdc', 'nlc'}``

            st = {'var', 'lin', 'qcn', 'nle', 'nli', 'qdc', 'nlc'};
        end

        function new_mm = copy(mm)
            % Duplicate the object.
            % ::
            %
            %   new_mm = mm.copy()
            %
            % Output:
            %   new_mm (mp.opt_model) : duplicate of original object
            %
            % Make a duplicate of the object, including duplicates of all
            % contained objects.

            %% delete old 'params' (cached parameters) fields
            %% to avoid clash with newer params() method
            fn = mm.get_set_types();
            for f = 1:length(fn)
                if isfield(mm.(fn{f}), 'params')
                    mm.(fn{f}) = rmfield(mm.(fn{f}), 'params');
                end
            end

            %% initialize copy
            new_mm = eval(class(mm));   %% create new object

            %% copy properties/fields
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(mm);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_mm = copy_prop(mm, new_mm, props{k});
            end
        end

        function s = to_struct(mm)
            % Convert object data *to* a struct.
            % ::
            %
            %   s = mm.to_struct()
            %
            % Converts the object data *to* a struct that can later be
            % converted back to an identical object using mp.struct2object.
            % Useful for saving the object data to a MAT-file in Octave.

            s = nested_struct_copy(struct(), mm);
            s.class_ = class(mm);
            st = mm.get_set_types();
            for k = 1:length(st)
                s.(st{k}) = s.(st{k}).to_struct();
            end
        end

        function mm = from_struct(mm, s)
            % Copy object data *from* a struct.
            % ::
            %
            %   mm.from_struct(s)
            %
            % Called by function mp.struct2object, after creating the object
            % to copy the object data *from* a struct. Useful for recreating
            % the object after loading struct data from a MAT-file in Octave.

            mm = nested_struct_copy(mm, s);
            st = mm.get_set_types();
            for k = 1:length(st)
                mm.(st{k}) = mp.struct2object(mm.(st{k}));
            end
        end

        function varargout = get_idx(mm, varargin)
            % Return ``idx`` struct for vars, constraints, costs.
            % ::
            %
            %   vv = mm.get_idx()
            %   [vv, ll] = mm.get_idx()
            %   [vv, ll, nne] = mm.get_idx()
            %   [vv, ll, nne, nni] = mm.get_idx()
            %   [vv, ll, nne, nni, qq] = mm.get_idx()
            %   [vv, ll, nne, nni, qq, nnc] = mm.get_idx()
            %   [vv, ll, nne, nni, qq, nnc, qqcn] = mm.get_idx()
            %
            % Returns the :attr:`idx <mp.set_manager.idx>` property of each
            % of the set manager objects, with the beginning and ending
            % index value and the number of elements for each named block.
            % The ``i1`` field is a struct with all of the  starting indices,
            % ``iN`` contains all the ending indices, and ``N`` contains all
            % the sizes. Each is a struct whose fields are the named blocks.
            %
            % Alternatively, you can specify the property name of the set type
            % directly as inputs.
            % ::
            %
            %   [idx1, idx2, ...] = mm.get_idx(set_type1, set_type2, ...)
            %
            % Inputs:
            %   set_type<n> (char array) : name of set type, valid options are:
            %
            %       - ``'var'`` - variables
            %       - ``'lin'`` - linear constraints
            %       - ``'qcn'`` - quadratic constraints
            %       - ``'nle'`` - nonlinear equality constraints
            %       - ``'nli'`` - nonlinear inequality constraints
            %       - ``'qdc'`` - quadratic costs
            %       - ``'nlc'`` - general nonlinear costs
            %
            % Outputs:
            %   vv (struct) : ``mm.var.idx``
            %   ll (struct) : ``mm.lin.idx``
            %   nne (struct) : ``mm.nle.idx``
            %   nni (struct) : ``mm.nli.idx``
            %   qq (struct) : ``mm.qdc.idx``
            %   nnc (struct) : ``mm.nlc.idx``
            %   qqcn (struct) : ``mm.qcn.idx``
            %   idx<n> (struct) : :attr:`idx <mp.set_manager.idx>` property
            %       of property indicated by respective input
            %
            % **Examples**::
            %
            %   [vv, ll, nne] = mm.get_idx();
            %   [vv, ll, qq] = mm.get_idx('var', 'lin', 'qdc');
            %   [ll, nne, nni] = mm.get_idx('lin', 'nle', 'nli');
            %
            % For a variable block named ``z`` we have ...
            %
            %   - ``vv.i1.z`` - starting index for ``z`` in vector ``x``
            %   - ``vv.iN.z`` - ending index for ``z`` in vector ``x``
            %   - ``vv.N.z`` - number of elements in ``z``
            %
            % To extract a ``z`` variable from ``x``::
            %
            %   z = x(vv.i1.z:vv.iN.z);
            %
            % To extract the multipliers on a linear constraint set
            % named ``foo``, where ``mu_l`` and ``mu_u`` are the full set of
            % linear constraint multipliers::
            %
            %   mu_l_foo = mu_l(ll.i1.foo:ll.iN.foo);
            %   mu_u_foo = mu_u(ll.i1.foo:ll.iN.foo);
            %
            % The number of nonlinear equality constraints in a set named
            % ``bar``::
            %
            %   nbar = nne.N.bar;
            %
            % *Note:* The following is preferable if you haven't already called
            % get_idx to get ``nne``.
            % ::
            %
            %   nbar = mm.nle.get_N('bar');
            %
            % If ``z``, ``foo``, and ``bar`` are indexed sets, then you can
            % replace them with something like ``z(i,j)``, ``foo(i,j,k)``
            % or ``bar(i)`` in the examples above.

            if nargin == 1
                varargout{1} = mm.var.idx;
                if nargout > 1
                    varargout{2} = mm.lin.idx;
                    if nargout > 2
                        varargout{3} = mm.nle.idx;
                        if nargout > 3
                            varargout{4} = mm.nli.idx;
                            if nargout > 4
                                varargout{5} = mm.qdc.idx;
                                if nargout > 5
                                    varargout{6} = mm.nlc.idx;
                                    if nargout > 6
                                        varargout{7} = mm.qcn.idx;
                                    end
                                end
                            end
                        end
                    end
                end
            else
                for k = nargout:-1:1
                    varargout{k} = mm.(varargin{k}).idx;
                end
            end
        end

        function rv = get_userdata(mm, name)
            % Used to retrieve values of user data.
            % ::
            %
            %   val = mm.get_userdata(name)
            %
            % Returns the value specified by the given name or an empty matrix
            % if userdata with ``name`` does not exist.
            %
            % Inputs:
            %   name (char array) : name of user data to return
            %
            % Outputs:
            %   val (arbitrary) : the arbitrary user data stored under the
            %       specified name
            %
            % This function allows the user to retrieve any arbitrary data that
            % was saved in the object for later use. Data for a given ``name``
            % is saved by assigning it to ``mm.userdata.(name)``.
            %
            % This can be useful, for example, when using a user function to add
            % variables or constraints, etc. Suppose some special indexing is
            % constructed when adding some variables or constraints. This
            % indexing data can be stored and used later to "unpack" the results
            % of the solved case.

            if isfield(mm.userdata, name)
                rv = mm.userdata.(name);
            else
                rv = [];
            end
        end

        function prob = problem_type(mm, recheck)
            % Return string identifying type of mathematical program.
            % ::
            %
            %   prob_type = mm.problem_type()
            %   prob_type = mm.problem_type(recheck)
            %
            % Returns a string identifying the type of mathematical program
            % represented by the current model, based on the variables, costs,
            % and constraints that have been added to the model. Used to
            % automatically select an appropriate solver.
            %
            % Linear and nonlinear equations are models with no costs, no
            % inequality constraints, and an equal number of continuous
            % variables and equality constraints. If the number of variables in
            % a nonlinear equation model is one more than the number of
            % constraints, it is a parameterized nonlinear equation.
            %
            % The output value is cached for future calls, but calling with a
            % true value for the optional ``recheck`` argument will force it to
            % recheck in case the problem type has changed due to modifying the
            % variables, constraints or costs in the model.
            %
            % Input:
            %   recheck (logical) : true to force a reevaluation instead of
            %       possibly returning a cached value
            %
            % Output:
            %     prob_type (char array) : problem type, one of the following
            %       strings:
            %
            %       - ``'LEQ'``   - linear equations
            %       - ``'NLEQ'``  - nonlinear equations
            %       - ``'PNE'``   - parameterized nonlinear equations
            %       - ``'LP'``    - linear program
            %       - ``'QP'``    - quadratic program
            %       - ``'NLP'``   - nonlinear program
            %       - ``'MILP'``  - mixed-integer linear program
            %       - ``'MIQP'``  - mixed-integer quadratic program
            %       - ``'MINLP'`` - mixed-integer nonlinear program

            if isempty(mm.prob_type) || nargin > 1 && recheck
                nleN = mm.nle.get_N();      %% nonlinear equalities
                nliN = mm.nli.get_N();      %% nonlinear inequalities
                nlcN = mm.nlc.get_N();      %% general nonlinear costs
                qdcN = mm.qdc.get_N();      %% quadratic costs
                qcnN = mm.qcn.get_N();      %% quadratic constraints
                linN = mm.lin.get_N();      %% linear constraints
                varN = mm.var.get_N();      %% variables
                if varN == 0
                    prob = '';
                elseif nlcN || qdcN         %% problem has costs
                    if nliN || nleN || nlcN %% nonlinear
                        prob = 'NLP';           %% nonlinear program
                    elseif qcnN                 %% quadratically constrained quadratic program
                        prob = 'QCQP';
                    else                    %% linear constraints, no general nonlinear costs
                        %% get quadratic cost coefficients
                        H = mm.qdc.params(mm.var);
                        if isempty(H) || ~any(any(H))
                            prob = 'LP';        %% linear program
                        else
                            prob = 'QP';        %% quadratic program
                        end
                    end
                else                    %% problem has no costs
                    if nliN
                        error('mp.opt_model.problem_type: invalid problem - nonlinear inequality constraints with no costs');
                    end
                    if nleN + qcnN + linN == varN || nleN + qcnN + linN == varN - 1   %% square (or almost) system
                        if linN > 0
                            %% get lower & upper bounds
                            [A, l, u] = mm.lin.params(mm.var);
                            if any(l ~= u)
                                error('mp.opt_model.problem_type: invalid problem - linear inequality constraints with no costs');
                            end
                        end
                        if nleN + qcnN + linN == varN  %% square system
                            if (nleN + qcnN) ~= 0
                                prob = 'NLEQ';      %% square nonlinear set of equations
                            else
                                prob = 'LEQ';       %% square linear set of equations
                            end
                        elseif nleN + linN + 1 == varN  %% square + 1 extra (parameterization) variable
                            if nleN
                                prob = 'PNE';       %% parameterized nonlinear set of equations
                            else
                                prob = 'PLEQ';      %% parameterized linear set of equations
                                error('mp.opt_model.problem_type: invalid problem - PNE not implemented for for linear constraints only');
                            end
                        else
                            error('mp.opt_model.problem_type: invalid problem - PNE must have num of vars = num of constraints + 1');
                        end
                    else
                        error('mp.opt_model.problem_type: invalid problem - non-square system with no costs');
                    end
                end
                if mm.is_mixed_integer() && ~strcmp(prob, 'NLEQ')
                    prob = ['MI' prob];
                end
                mm.prob_type = prob;    %% cache it
            else
                prob = mm.prob_type;    %% return cached type
            end
        end

        function TorF = is_mixed_integer(mm)
            % Return true if model is mixed integer, false otherwise.
            % ::
            %
            %   TorF = mm.is_mixed_integer()
            %
            % Outputs:
            %   TorF (logical): true or false, indicating whether any of the
            %       variables are binary or integer

            TorF = 0;
            if mm.var.get_N()
                for k = 1:length(mm.var.order)
                    t = mm.var.data.vt.(mm.var.order(k).name);
                    if iscell(t)
                        for j = 1:length(t(:))
                            if any(t{j} ~= 'C')
                                TorF = 1;
                                break;
                            end
                        end
                    else
                        if any(t ~= 'C')
                            TorF = 1;
                            break;
                        end
                    end
                end
            end
        end

        function TorF = is_solved(mm)
            % Return true if model has been solved.
            % ::
            %
            %   TorF = mm.is_solved()
            %
            % Outputs:
            %   TorF (logical): true or false, indicating whether the model
            %       solution is available

            TorF = ~isempty(mm.soln.eflag);
        end

        function [x, f, eflag, output, lambda] = solve(mm, opt)
            % Solve the model.
            % ::
            %
            %   x = mm.solve()
            %   [x, f] = mm.solve()
            %   [x, f, exitflag] = mm.solve()
            %   [x, f, exitflag, output] = mm.solve()
            %   [x, f, exitflag, output, jac] = mm.solve()      (LEQ/NLEQ problems)
            %   [x, f, exitflag, output, lambda] = mm.solve()   (other problem types)
            %   [x ...] = mm.solve(opt)
            %
            % Solves the model using one of the following, depending on the
            % problem type: mp_linsolve, nleqs_master, pnes_master,
            % qps_master, qcqps_master, miqps_master, nlps_master.
            %
            % Inputs:
            %   opt (struct) : optional options struct with the following
            %       fields, all of which are also optional *(default values
            %       shown in parentheses)*
            %
            %       - ``alg`` (``'DEFAULT'``) - determines which solver to
            %         use, list of relevant problem types are listed in parens
            %         next to each
            %
            %           - ``'DEFAULT'`` - automatic, depending on problem type,
            %             uses the the first available of:
            %
            %             ===== ================================================
            %             Type  Solver Precedence
            %             ===== ================================================
            %             LP    Gurobi, CPLEX, MOSEK, linprog *(if MATLAB)*,
            %                   HIGHS, GLPK, BPMPD, MIPS
            %             QP    Gurobi, CPLEX, MOSEK, quadprog *(if MATLAB)*,
            %                   HIGHS, BPMPD, MIPS
            %             QCQP  IPOPT, Artelys Knitro, FMINCON, MIPS
            %             MILP  Gurobi, CPLEX, MOSEK, Opt Tbx (intlingprog),
            %                   HIGHS, GLPK
            %             MIQP  Gurobi, CPLEX, MOSEK
            %             NLP   MIPS
            %             MINLP Artelys Knitro *(not yet implemented)*
            %             LEQ   built-in backslash operator
            %             NLEQ  Newton's method
            %             PNE   predictor/corrector method
            %             ===== ================================================
            %
            %           - ``'BPMPD'``   - *(LP, QP)* BPMPD_MEX
            %           - ``'CLP'``     - *(LP, QP)* CLP
            %           - ``'CPLEX'``   - *(LP, QP, MILP, MIQP)* CPLEX
            %           - ``'FD'``      - *(NLEQ)* fast-decoupled Newon's method
            %           - ``'FMINCON'`` - *(QCQP, NLP)* FMINCON, MATLAB
            %             Optimization Toolbox
            %           - ``'FSOLVE'``  - *(NLEQ)* FSOLVE, MATLAB Optimization
            %             Toolbox
            %           - ``'GLPK'``    - *(LP, MILP)* GLPK
            %           - ``'GS'``      - *(NLEQ)* Gauss-Seidel
            %           - ``'GUROBI'``  - *(LP, QP, QCQP, MILP, MIQP)* Gurobi
            %           - ``'HIGHS'``   - *(LP, QP, MILP)* HiGHS,
            %             https://highs.dev
            %           - ``'IPOPT'``   - *(LP, QP, QCQP, NLP)* IPOPT, requires
            %             MEX interface to IPOPT solver
            %             https://github.com/coin-or/Ipopt
            %           - ``'KNITRO'``  - *(QCQP, NLP, MINLP)* Artelys Knitro,
            %             requires Artelys Knitro solver
            %             https://www.artelys.com/solvers/knitro/
            %           - ``'MIPS'``    - *(LP, QP, QCQP, NLP)* MIPS, MATPOWER
            %             Interior Point Solver, pure MATLAB implementation of
            %             a primal-dual interior point method, if
            %             ``mips_opt.step_control = 1`` it uses MIPS-sc, a
            %             step controlled variant of MIPS
            %           - ``'MOSEK'``   - *(LP, QP, MILP, MIQP)* MOSEK
            %           - ``'NEWTON'``  - *(NLEQ)* Newton's method
            %           - ``'OSQP'``    - *(LP, QP)* OSQP, https://osqp.org
            %           - ``'OT'``      - *(LP, QP, MILP)* MATLAB Optimization
            %             Toolbox, linprog, quadprog or intlinprog
            %       - ``verbose`` (0) - controls level of progress output
            %         displayed:
            %
            %           - 0 = no progress output
            %           - 1 = some progress output
            %           - 2 = verbose progress output
            %       - ``bp_opt``      - options vector for :func:`bp` (BPMPD)
            %       - ``clp_opt``     - options vector for CLP
            %       - ``cplex_opt``   - options struct for CPLEX
            %       - ``fd_opt``      - options struct for fast-decoupled
            %         Newton's method, nleqs_fd_newton
            %       - ``fmincon_opt`` - options struct for :func:`fmincon`
            %       - ``fsolve_opt``  - options struct for :func:`fsolve`
            %       - ``glpk_opt``    - options struct for GLPK
            %       - ``grb_opt``     - options struct for Gurobi
            %       - ``gs_opt``      - options struct for Gauss-Seidel method,
            %         nleqs_gauss_seidel
            %       - ``highs_opt``   - options struct for HiGHS
            %       - ``intlinprog_opt`` - options struct for :func:`intlinprog`
            %       - ``ipopt_opt``   - options struct for IPOPT
            %       - ``knitro_opt``  - options struct for Artelys Knitro
            %       - ``leq_opt``     - options struct for mplinsolve, with
            %         optional fields ``solver`` and ``opt`` corresponding to
            %         respective mplinsolve args, and ``thresh`` specifying a
            %         threshold on the absolute value of any element ``x``,
            %         above which ``exitflag`` will be set to 0
            %       - ``linprog_opt`` - options struct for :func:`linprog`
            %       - ``mips_opt``    - options struct for mips
            %       - ``mosek_opt``   - options struct for MOSEK
            %       - ``newton_opt``  - options struct for Newton method,
            %         nleqs_newton
            %       - ``osqp_opt``    - options struct for OSQP
            %       - ``quadprog_opt`` - options struct for :func:`quadprog`
            %       - ``parse_soln`` (0) - flag that specifies whether or not
            %         to call the parse_soln() method and place the return values
            %         in ``mm.soln``.
            %       - ``price_stage_warn_tol`` (1e-7) - tolerance on 
            %         objective fcn value and primal variable relative match
            %         required to avoid mismatch warning message if mixed
            %         integer price computation stage is not skipped
            %       - ``relax_integer`` (0) - relax integer constraints, if true
            %       - ``skip_prices`` (0) - flag that specifies whether or not
            %         to skip the price computation stage for mixed integer
            %         problems, in which the problem is re-solved for only the
            %         continuous variables, with all others being constrained
            %         to their solved values
            %       - ``x0`` (empty)  - optional initial value of ``x``,
            %         overrides value stored in model *(ignored by some solvers)*
            %
            % Outputs:
            %     x (double) : solution vector
            %     f (double) : final (objective) function value
            %     exitflag (integer) : exit flag
            %
            %       - 1 = converged
            %       - 0 or negative values = solver specific failure codes
            %     output (struct) : output struct with the following fields:
            %
            %       - ``alg`` - algorithm code of solver used
            %       - ``et``  - elapsed time (sec)
            %       - *(others)* - solver specific fields
            %     jac (double) : final Jacobian matrix *(if available, for
            %       LEQ/NLEQ problems)*
            %     lambda (struct) : *(for all problem types except LEQ, NLEQ,
            %       and PNE)* Langrange and Kuhn-Tucker multipliers on the
            %       constraints, struct with fields:
            %
            %       - ``eqnonlin`` - nonlinear equality constraints
            %       - ``ineqnonlin`` - nonlinear inequality constraints
            %       - ``mu_l`` - lower (left-hand) limit on linear constraints
            %       - ``mu_u`` - upper (right-hand) limit on linear constraints
            %       - ``lower`` - lower bound on optimization variables
            %       - ``upper`` - upper bound on optimization variables
            %
            % See also mp_linsolve, nleqs_master, pnes_master, qps_master,
            % qcqps_master, miqps_master, nlps_master.

            t0 = tic;       %% start timer
            if nargin < 2
                opt = struct();
            end
            % opt.parse_soln = 1;

            %% call appropriate solver
            pt = mm.problem_type();
            switch pt
                case 'LEQ'          %% LEQ   - linear equations
                    if isfield(opt, 'leq_opt')
                        if isfield(opt.leq_opt, 'solver')
                            leq_solver = opt.leq_opt.solver;
                        else
                            leq_solver = '';
                        end
                        if isfield(opt.leq_opt, 'opt')
                            leq_opt = opt.leq_opt.opt;
                        else
                            leq_opt = struct();
                        end
                        if isfield(opt.leq_opt, 'thresh')
                            leq_thresh = opt.leq_opt.thresh;
                        else
                            leq_thresh = 0;
                        end
                    else
                        leq_solver = '';
                        leq_opt = struct();
                        leq_thresh = 0;
                    end

                    [A, b] = mm.lin.params(mm.var);
                    if leq_thresh           %% check for failure
                        %% set up to trap non-singular matrix warnings
                        [lastmsg, lastid] = lastwarn;
                        lastwarn('');

                        x = mplinsolve(A, b, leq_solver, leq_opt);

                        [msg, id] = lastwarn;
                        %% Octave is not consistent in assigning proper warning id,
                        %% so we just check for presence of *any* warning
                        if ~isempty(msg) || max(abs(x)) > leq_thresh
                            eflag = 0;
                        else
                            eflag = 1;
                        end
                    else                    %% no failure check
                        x = mplinsolve(A, b, leq_solver, leq_opt);
                        eflag = 1;
                    end
                    f = A*x - b;
                    output = struct('alg', leq_solver);
                    lambda = A;     %% jac
                case {'NLEQ', 'PNE'}    %% NLEQ, PNE - nonlinear equations
                    if isfield(opt, 'x0')
                        x0 = opt.x0;
                    else
                        x0 = mm.var.params();
                    end

                    fcn = @(x)nleq_fcn_(mm, x);
                    switch pt
                        case 'NLEQ' %% NLEQ - nonlinear equation
                            [x, f, eflag, output, lambda] = nleqs_master(fcn, x0, opt);
                        case 'PNE'  %% PNE - parameterized nonlinear equation
                            [x, f, eflag, output, lambda] = pnes_master(fcn, x0, opt);
                    end
                case {'MINLP', 'NLP'}
                    mixed_integer = strcmp(pt(1:2), 'MI') && ...
                        (~isfield(opt, 'relax_integer') || ~opt.relax_integer);
                    if mixed_integer    %% MINLP - mixed integer non-linear program
                        error('mp.opt_model.solve: not yet implemented for ''MINLP'' problems.')
                    else                %% NLP   - nonlinear program
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax] = mm.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end

                        %% run solver
                        [A, l, u] = mm.lin.params(mm.var);
                        f_fcn = @(x)nlp_costfcn(mm, x);
                        gh_fcn = @(x)nlp_consfcn(mm, x);
                        hess_fcn = @(x, lambda, cost_mult)nlp_hessfcn(mm, x, lambda, cost_mult);
                        [x, f, eflag, output, lambda] = ...
                            nlps_master(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
                    end
                otherwise
                    %% get parameters
                    [HH, CC, C0] = mm.qdc.params(mm.var);
                    [Q, B, ll, uu] = mm.qcn.params(mm.var);
                    [A, l, u] = mm.lin.params(mm.var);
                    mixed_integer = strcmp(pt(1:2), 'MI') && ...
                        (~isfield(opt, 'relax_integer') || ~opt.relax_integer);

                    if mixed_integer
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax, vtype] = mm.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end

                        %% run solver
                        if isempty(Q)          %% MILP, MIQP - mixed integer linear/quadratic program
                            [x, f, eflag, output, lambda] = ...
                                miqps_master(HH, CC, A, l, u, xmin, xmax, x0, vtype, opt);
                        else                   %% MIQCQP - mixed integer quadratically constrained quadratic program
                            % To be implemented ...
                            % [x, f, eflag, output, lambda] = ...
                            %    miqcqps_master(HH, CC, Q, B, k, ll, uu, A, l, u, xmin, xmax, x0, vtype, opt);
                        end
                    else                %% LP, QP - linear/quadratic program
                        %% optimization vars, bounds, types
                        [x0, xmin, xmax] = mm.var.params();
                        if isfield(opt, 'x0')
                            x0 = opt.x0;
                        end

                        %% run solver
                        if isempty(Q)          %% LP, QP - linear/quadratic program
                            [x, f, eflag, output, lambda] = ...
                                qps_master(HH, CC, A, l, u, xmin, xmax, x0, opt);
                        else                   %% QCQP - quadratically constrained quadratic program
                            [x, f, eflag, output, lambda] = ...
                                qcqps_master(HH, CC, Q, B, ll, uu, A, l, u, xmin, xmax, x0, opt);
                        end
                    end
                    f = f + C0;
            end

            %% store solution
            mm.soln.eflag = eflag;
            mm.soln.x = x;
            mm.soln.f = f;
            mm.soln.output = output;
            if isstruct(lambda)
                mm.soln.lambda = lambda;
            else
                mm.soln.jac = lambda;
            end

            %% parse solution
            if isfield(opt, 'parse_soln') && opt.parse_soln
                mm.parse_soln(true);
            end
            mm.soln.output.et = toc(t0);    %% stop timer
        end

        function TorF = has_parsed_soln(mm)
            % Return true if model has a parsed solution.
            % ::
            %
            %   TorF = mm.has_parsed_soln()

            TorF = mm.var.has_parsed_soln();
        end

        function ps = parse_soln(mm, stash)
            % Parse solution vector and shadow prices by all named sets.
            % ::
            %
            %   ps = mm.parse_soln()
            %   mm.parse_soln(stash)
            %
            % For a solved model, parse_soln returns a struct of parsed solution
            % vector and shadow price values for each named set of variables
            % and constraints.
            %
            % Input:
            %   stash (logical) : if true, individual parsed solutions are
            %       stored in the :attr:`soln` property of the respective
            %       set type object.
            %
            % Output:
            %   ps (struct) : struct of parsed solution values with the
            %       following fields, where each of the terminal elements is
            %       a struct with fields corresponding to the respective named
            %       sets.
            %
            %       - ``var`` - variable solution struct with fields:
            %
            %           - ``val`` - values of variables
            %           - ``mu_l`` - shadow prices on variable lower bounds
            %           - ``mu_u`` - shadow prices on variable upper bounds
            %       - ``lin`` - linear constraint solution struct with fields:
            %
            %           - ``mu_l`` - shadow prices on constraint lower bounds
            %           - ``mu_u`` - shadow prices on constraint upper bounds
            %       - ``qcn`` - quadratic constraint solution struct with fields:
            %
            %           - ``mu_l`` - shadow prices on constraint lower bounds
            %           - ``mu_u`` - shadow prices on constraint upper bounds
            %       - ``nle`` - nonlinear equality constraint solution struct
            %         with fields:
            %
            %           - ``lam`` - shadow prices on constraints
            %       - ``nli`` - nonlinear inequality constraint solution struct
            %         with fields:
            %
            %           - ``mu`` - shadow prices on constraints
            %
            % The value of each element in the returned struct can be obtained
            % via the :meth:`get_soln` method of the appropriate set type
            % object as well, but using parse_soln is generally more efficient
            % if a complete set of values is needed.

            if nargin < 2
                stash = false;
            end

            if ~mm.is_solved()
                error('mp.opt_model.parse_soln: model not solved');
            end

            %% var
            ps = struct('var', mm.var.parse_soln(mm.soln, stash));

            %% lin
            ps_lin = mm.lin.parse_soln(mm.soln, stash);
            if ~isempty(ps_lin)
                ps.lin = ps_lin;
            end

            %% qcn
            ps_qcn = mm.qcn.parse_soln(mm.soln, stash);
            if ~isempty(ps_qcn)
                ps.qcn = ps_qcn;
            end

            %% nle
            ps_nle = mm.nle.parse_soln(mm.soln, true, stash);
            if ~isempty(ps_nle)
                ps.nle = ps_nle;
            end

            %% nli
            ps_nli = mm.nli.parse_soln(mm.soln, false, stash);
            if ~isempty(ps_nli)
                ps.nli = ps_nli;
            end

            %%-----  DEPRECATED  -----
            %% if requested, stash the result directly in mm.soln
            %% (they are already stashed in the soln property of each set type)
            % if stash
            %     mm.soln = nested_struct_copy(mm.soln, ps, struct('copy_mode', '='));
            % end
        end

        function display(mm, varargin)
            % Display the object.
            %
            % ::
            %
            %   mm
            %
            % Called when semicolon is omitted at the command-line. Displays the details
            % of the variables, constraints, costs included in the model.

            if nargin < 2
                more_set_types = {};
            else
                more_set_types = varargin{1};
            end

            mm.display_header();

            %% display details of each set type
            set_types = mm.get_set_types();
            set_types = horzcat(set_types, more_set_types);
            % set_types = {'var', 'nle', 'nli', 'lin', 'qcn', 'qdc', 'nlc', more_set_types{:}};
            fprintf('\n');
            for k = 1:length(set_types)
                mm.(set_types{k}).display(set_types{k});
            end

            mm.display_footer();
        end

        function mm = display_header(mm)
            % Display header, before each set type.
            % ::
            %
            %   mm.display_header()
            %
            % Called automatically by display, *before* displaying each
            % set type.

            fprintf('CLASS : %s\n', class(mm));
        end

        function mm = display_footer(mm)
            % Display footer, after each set type.
            % ::
            %
            %   mm.display_footer()
            %
            % Called automatically by display, *after* displaying each
            % set type. By default, displays information about user data.

            %% user data
            fields = fieldnames(mm.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                               size       class\n');
                fprintf(' ------------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = mm.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-31s %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            else
                fprintf('USER DATA                   :  <none>\n');
            end
        end

        function mm = display_soln(mm, varargin)
            % Display solution values.
            % ::
            %
            %   mm.display_soln()
            %   mm.display_soln(set_type)
            %   mm.display_soln(set_type, name)
            %   mm.display_soln(set_type, name, idx)
            %   mm.display_soln(fid)
            %   mm.display_soln(fid, set_type)
            %   mm.display_soln(fid, set_type, name)
            %   mm.display_soln(fid, set_type, name, idx)
            %
            % Displays the model solution, including values, bounds and shadow
            % prices for variables, linear and quadratic constraints, values
            % and shadow prices for nonlinear constraints, and individual cost
            % components.
            %
            % Results are displayed for all set types or those specified by
            % ``set_type`` and for each named/indexed set or a specified
            % ``name``/``idx``.
            %
            % Inputs:
            %   set_type (char array or cell array) : one of the following,
            %       specifying the type of set:
            %
            %       - ``'var'`` - variables
            %       - ``'lin'`` - linear constraints
            %       - ``'qcn'`` - quadratic constraints
            %       - ``'nle'`` - nonlinear equality constraints
            %       - ``'nli'`` - nonlinear inequality constraints
            %       - ``'nlc'`` - nonlinear costs
            %       - ``'qdc'`` - quadratic costs
            %
            %       *or*
            %
            %       - a cell array of one or more of the above
            %
            %       *or*
            %
            %       - ``''`` or ``'all'`` - indicating to display all
            %   name (char array) : *(optional)* the name of the set
            %   idx (cell array) : *(optional)* the indices of the set
            %
            % Examples::
            %
            %   mm.display_soln('var');
            %   mm.display_soln({'nle', 'nli'});
            %   mm.display_soln('var', 'P');
            %   mm.display_soln('lin', 'lin_con_1');
            %   mm.display_soln('nle', 'nle_con_b', {2,3});
            %
            % See also parse_soln.

            %% input arg handling
            if nargin < 2 || ischar(varargin{1})
                fid = 1;
                args = varargin;
            else
                fid = varargin{1};
                args = varargin(2:end);
            end
            nargs = length(args);

            set_type = 'all';
            name = [];
            idx = [];
            if nargs >= 1
                set_type = args{1};
                if nargs >= 2
                    name = args{2};
                    if nargs >= 3
                        idx = args{3};
                    end
                end
            end

            %% print header
            if mm.is_solved()
                if strcmp(set_type, 'all')
                    set_types = mm.get_set_types(); %% all set types
                elseif ~iscell(set_type)
                    set_types = {set_type}; %% make set_type a cell array of char arrays
                else
                    set_types = set_type;
                end

                for ss = 1:length(set_types)
                    st = set_types{ss};
                    sm = mm.(st);

                    switch st
                    case 'var'
                        sm.display_soln(mm.soln, fid, args{2:end});
                    case 'nle'
                        sm.display_soln(mm.var, mm.soln, 1, fid, args{2:end});
                    case 'nli'
                        sm.display_soln(mm.var, mm.soln, 0, fid, args{2:end});
                    otherwise
                        sm.display_soln(mm.var, mm.soln, fid, args{2:end});
                    end
                end             %% loop over set types
            else
                fprintf(fid, 'Not a solved model.\n');
            end
        end
    end     %% methods
end         %% classdef

function d = copy_prop(s, d, prop)
    %
    if isa(s.(prop), 'mp.sm_quad_cost_legacy')
        d.(prop) = s.(prop).copy('mp.sm_quad_cost');
    elseif isa(s.(prop), 'mp.set_manager')
        d.(prop) = s.(prop).copy();
    elseif isa(d.(prop), 'mp.set_manager')
        d.(prop) = nested_struct_copy( ...
            d.(prop), s.(prop));
    else
        d.(prop) = s.(prop);
    end
end

%% system of nonlinear and linear equations
function [f, J] = nleq_fcn_(mm, x)
    %
    flin = []; Jlin = [];
    fqcn = []; Jqcn = [];
    fnln = []; Jnln = [];
    if mm.lin.get_N()
        [flin, ~, Jlin] = mm.lin.eval(mm.var, x);
    end
    if nargout > 1
        if mm.qcn.get_N()
            [fqcn, Jqcn] = mm.qcn.eval(mm.var, x);
        end
        if mm.nle.get_N()
            [fnln, Jnln] = mm.nle.eval(mm.var, x);
        end
        J = [Jnln; Jqcn; Jlin];
    else
        if mm.qcn.get_N()
            fqcn = mm.qcn.eval(mm.var, x);
        end
        if mm.nle.get_N()
            fnln = mm.nle.eval(mm.var, x);
        end
    end
    f = [fnln; fqcn; flin];
end
