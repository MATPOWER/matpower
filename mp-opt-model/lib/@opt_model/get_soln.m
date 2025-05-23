function varargout = get_soln(om, set_type, varargin)
% get_soln - Fetch solution values for specific named/indexed sets.
% ::
%
%   VALS = OM.GET_SOLN(SET_TYPE, NAME)
%   VALS = OM.GET_SOLN(SET_TYPE, NAME, IDX)
%   VALS = OM.GET_SOLN(SET_TYPE, TAGS, NAME)
%   VALS = OM.GET_SOLN(SET_TYPE, TAGS, NAME, IDX)
%
%   Returns named/indexed results for a solved model, evaluated at
%   the solution found.
%
%   Inputs:
%       SET_TYPE - one of the following, specifying the type of set:
%           'var' - variables
%           'lin' - linear constraints
%           'qcn' - quadratic constraints
%           'nle' - nonlinear equality constraints
%           'nli' - nonlinear inequality constraints
%           'nlc' - nonlinear costs
%           'qdc' - quadratic costs
%       TAGS - char array or cell array of char arrays specifying the
%           desired output(s). Valid tags vary by SET_TYPE as follows:
%           'var' - default is {'x', 'mu_l', 'mu_u'}
%               'x' - value of solution variable
%               'mu_l' - shadow price on variable lower bound
%               'mu_u' - shadow price on variable upper bound
%           'lin' - default is {'g', 'mu_l', 'mu_u'}
%               'g' - 1 x 2 cell array of upper and lower constraint
%                   values, {A*x - u, l - A*x}
%               'Ax_u' - upper constraint value, A*x - u
%               'l_Ax' - lower constraint value, l - A*x
%               'mu_l' - shadow price on constraint lower bound
%               'mu_u' - shadow price on constraint upper bound
%           'qcn' - default is {'g', 'mu_l', 'mu_u'}
%               'g' - 1 x 2 cell array of upper and lower constraint
%                   values, {g(x) - u, l - g(x)}
%               'g_u' - upper constraint value, g(x) - u
%               'l_g' - lower constraint value, l - g(x)
%               'mu_l' - shadow price on constraint lower bound
%               'mu_u' - shadow price on constraint upper bound
%           'nle' - default is {'g', 'lam', 'dg'}
%               'g' - constraint value g(x)
%               'lam' - shadow price on constraint
%               'dg' - Jacobian of constraint
%           'nli' - default is {'h', 'mu', 'dh'}
%               'h' - constraint value h(x)
%               'mu' - shadow price on constraint
%               'dh' - Jacobian of constraint
%           'nlc' and 'qdc' - default is {'f', 'df', 'd2f'}
%               'f' - cost function value f(x) (for 'qdc' can return a vector)
%               'df' - gradient of cost function
%               'd2f' - Hessian of cost function
%       NAME - char array specifying the name of the set
%       IDX  - cell array specifying the indices of the set
%
%   Outputs:
%       Variable number of outputs corresponding to TAGS input. If TAGS
%       is empty or not specified, the calling context will define the
%       number of outputs, returned in order of default tags for the
%       specified SET_TYPE.
%
%   Examples:
%       [P, muPmin, muPmax] = om.get_soln('var', 'P');
%       [mu_u, mu_l] = om.get_soln('lin', {'mu_u', 'mu_l'}, 'lin_con_1');
%       dg_b_2_3 = om.get_soln('nle', 'dg', 'nle_con_b', {2,3});
%
%   For a complete set of solution vector values and shadow prices, using
%   the PARSE_SOLN method may be more efficient.
%
% See also parse_soln.

%   MP-Opt-Model
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

switch set_type
    case 'var'
        [varargout{1:nargout}] = om.var.get_soln(om.soln, varargin{:});
    case 'nle'
        [varargout{1:nargout}] = om.nle.get_soln(om.var, om.soln, true, varargin{:});
    case 'nli'
        [varargout{1:nargout}] = om.nli.get_soln(om.var, om.soln, false, varargin{:});
    case {'nlc', 'qdc'}
        [varargout{1:nargout}] = om.(set_type).get_soln(om.var, om.soln, varargin{:});
    case 'lin'
        if strcmp(om.problem_type(), 'LEQ')
            %% tag must be 'f'
            if length(varargin) >= 2 && ischar(varargin{2})
                varargin{1} = {'f'};
            else
                varargin = {'f', varargin{:}};
            end
        end
        [varargout{1:nargout}] = om.lin.get_soln(om.var, om.soln, varargin{:});
    case 'qcn'
        if strcmp(om.problem_type(), 'NLEQ')
            %% tag must be 'f'
            if length(varargin) >= 2 && ischar(varargin{2})
                varargin{1} = {'f'};
            else
                varargin = {'f', varargin{:}};
            end
        end
        [varargout{1:nargout}] = om.lin.get_soln(om.var, om.soln, varargin{:});
end
