function varargout = get_soln(om, set_type, tags, name, idx)
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

%% input arg handling
if nargin == 3              %% om.get_soln(set_type, name)
    idx = [];
    name = tags;
    tags = {};
elseif nargin == 4
    if ischar(name)         %% om.get_soln(set_type, tags, name)
        idx = [];
    else                    %% om.get_soln(set_type, name, idx)
        idx = name;
        name = tags;
        tags = {};
    end
end

%% set up tags for default outputs
if isempty(tags)
    switch set_type
        case 'var'
            tags = {'x', 'mu_l', 'mu_u'};
        case 'lin'
            if strcmp(om.problem_type(), 'LEQ')
                tags = {'f'};   %% 'Ax_u', 'l_Ax' are also options
            else
                tags = {'g', 'mu_l', 'mu_u'};   %% 'Ax_u', 'l_Ax' are also options
            end
        case 'nle'
            tags = {'g', 'lam', 'dg'};
        case 'nli'
            tags = {'h', 'mu', 'dh'};
        case 'nlc'
            tags = {'f', 'df', 'd2f'};
        case 'qdc'
            tags = {'f', 'df', 'd2f'};
    end
elseif ~iscell(tags)
    tags = { tags };
end

%% set up indexing
om_ff = om.(set_type);
if isempty(idx)         %% simple named set
    N = om_ff.idx.N.(name);
    i1 = om_ff.idx.i1.(name);           %% starting row index
    iN = om_ff.idx.iN.(name);           %% ending row index
else                    %% indexed named set
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields, updating only
    %% the subscripts before use
    sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
    sn = sc; sn(2).type = '()';                         %% num array field
    N = subsref(om_ff.idx.N, sn);
    i1 = subsref(om_ff.idx.i1, sn);     %% starting row index
    iN = subsref(om_ff.idx.iN, sn);     %% ending row index
end

%% get outputs
varargout = cell(1, nargout);
s = om.soln;
if N && ~isempty(s.eflag)
    switch set_type
        case 'var'
            for k = 1:nargout
                switch tags{k}
                    case 'x'
                        varargout{k} = s.x(i1:iN);
                    case 'mu_l'
                        varargout{k} = s.lambda.lower(i1:iN);
                    case 'mu_u'
                        varargout{k} = s.lambda.upper(i1:iN);
                    otherwise
                        error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                end
            end
        case 'lin'
            if strcmp(om.problem_type(), 'LEQ') %% tag must be 'f'
                varargout{1} = s.f(i1:iN);
            else
                if any(ismember({'g', 'Ax_u', 'l_Ax'}, tags(1:nargout)))
                    g = cell(1,2);
                    [g{:}] = om.eval_lin_constraint(s.x, name, idx);
                end
                for k = 1:nargout
                    switch tags{k}
                        case 'g'
                            varargout{k} = g;
                        case 'Ax_u'
                            varargout{k} = g{1};
                        case 'l_Ax'
                            varargout{k} = g{2};
                        case 'mu_l'
                            varargout{k} = s.lambda.mu_l(i1:iN);
                        case 'mu_u'
                            varargout{k} = s.lambda.mu_u(i1:iN);
                        otherwise
                            error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                    end
                end
            end
        case 'nle'
            if ismember('dg', tags(1:nargout))
                [g, dg] = om.eval_nln_constraint(s.x, 1, name, idx);
            elseif ismember('g', tags(1:nargout))
                g = om.eval_nln_constraint(s.x, 1, name, idx);
            end
            for k = 1:nargout
                switch tags{k}
                    case 'g'
                        varargout{k} = g;
                    case 'dg'
                        varargout{k} = dg;
                    case 'lam'
                        varargout{k} = s.lambda.eqnonlin(i1:iN);
                    otherwise
                        error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                end
            end
        case 'nli'
            if ismember('dh', tags(1:nargout))
                [h, dh] = om.eval_nln_constraint(s.x, 0, name, idx);
            elseif ismember('h', tags(1:nargout))
                h = om.eval_nln_constraint(s.x, 0, name, idx);
            end
            for k = 1:nargout
                switch tags{k}
                    case 'h'
                        varargout{k} = h;
                    case 'dh'
                        varargout{k} = dh;
                    case 'mu'
                        varargout{k} = s.lambda.ineqnonlin(i1:iN);
                    otherwise
                        error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                end
            end
        case 'nlc'
            if ismember('d2f', tags(1:nargout))
                [f, df, d2f] = om.eval_nln_cost(s.x, name, idx);
            elseif ismember('df', tags(1:nargout))
                [f, df] = om.eval_nln_cost(s.x, name, idx);
            else
                f = om.eval_nln_cost(s.x, name, idx);
            end
            for k = 1:nargout
                switch tags{k}
                    case 'f'
                        varargout{k} = f;
                    case 'df'
                        varargout{k} = df;
                    case 'd2f'
                        varargout{k} = d2f;
                    otherwise
                        error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                end
            end
        case 'qdc'
            if ismember('d2f', tags(1:nargout))
                [f, df, d2f] = om.eval_quad_cost(s.x, name, idx);
            elseif ismember('df', tags(1:nargout))
                [f, df] = om.eval_quad_cost(s.x, name, idx);
            else
                f = om.eval_quad_cost(s.x, name, idx);
            end
            for k = 1:nargout
                switch tags{k}
                    case 'f'
                        varargout{k} = f;
                    case 'df'
                        varargout{k} = df;
                    case 'd2f'
                        varargout{k} = d2f;
                    otherwise
                        error('opt_model.get_soln: unknown tag ''%s''', tags{k});
                end
            end
        otherwise
            error('opt_model.get_soln: unknown set_type ''%s''', set_type);
    end     %% switch set_type
    end
end     %% if N
