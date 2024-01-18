classdef cost_table < mp_table_subclass
% mp.cost_table - Table for (polynomial and piecewise linear) cost parameters.
% ::
%
%   T = cost_table(poly_n, poly_coef, pwl_n, pwl_qty, pwl_cost);
%
% .. important::
%
%   Since the dot syntax ``T.<var_name>`` is used to access table variables,
%   you must use a functional syntax ``<method>(T,...)``, as opposed to
%   the object-oriented ``T.<method>(...)``, to call standard
%   mp.cost_table methods.
%
% Standard table subscripting syntax is not available within methods of
% this class (references built-in :func:`subsref` and :func:`subsasgn`
% rather than the versions overridden by the table class). For this reason,
% some method implementations are delegated to static methods in
% mp.cost_table_utils where that syntax is available, making the code more
% readable.
%
% mp.cost_table Methods:
%   * cost_table - construct object
%   * poly_params - create struct of polynomial parameters from mp.cost_table
%   * pwl_params - create struct of piecewise linear parameters from mp.cost_table
%   * max_pwl_cost - get maximum cost component used to specify pwl costs
%
% An mp.cost_table has the following columns:
%
%   =============  =========  =====================================
%   Name           Type       Description
%   =============  =========  =====================================
%   ``poly_n``     *integer*  :math:`n_\mathrm{poly}`, number of coefficients
%                             in polynomial cost curve,
%                             :math:`f_\mathrm{poly}(x) = c_0 + c_1 x... + c_N x^N`,
%                             where :math:`n_\mathrm{poly} = N + 1`
%   ``poly_coef``  *double*   matrix of coefficients :math:`c_j`, of polynomial
%                             cost :math:`f_\mathrm{poly}(x)`, where :math:`c_j`
%                             is found in column :math:`j+1`
%   ``pwl_n``      *double*   :math:`n_\mathrm{pwl}`, number of data points
%                             :math:`(x_1, f_1), (x_2, f_2), ..., (x_N, f_N)`
%                             defining a piecewise linear cost curve,
%                             :math:`f_\mathrm{pwl}(x)` where
%                             :math:`N = n_\mathrm{pwl}`
%   ``pwl_qty``    *double*   matrix of *quantiy* coordinates :math:`x_j` for
%                             piecwise linear cost :math:`f_\mathrm{pwl}(x)`,
%                             where :math:`x_j` is found in column :math:`j`
%   ``pwl_cost``   *double*   matrix of *cost* coordinates :math:`f_j` for
%                             piecwise linear cost :math:`f_\mathrm{pwl}(x)`,
%                             where :math:`f_j` is found in column :math:`j`
%   =============  =========  =====================================
%
% See also mp.cost_table_utils, mp_table_subclass.

%   MATPOWER
%   Copyright (c) 2023-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = cost_table(varargin)
            % ::
            %
            %   T = cost_table()
            %   T = cost_table(poly_n, poly_coef, pwl_n, pwl_qty, pwl_cost)
            %
            % *For descriptions of the inputs, see the corresponding column in
            % the class documentation above.*
            %
            % Inputs:
            %   poly_n (col vector of integers)
            %   poly_coef (matrix of doubles)
            %   pwl_n (col vector of integers)
            %   pwl_qty (matrix of doubles)
            %   pwl_cost (matrix of doubles)
            %
            % Outputs:
            %   T (mp.cost_table) : the cost table object

            if nargin == 0
                args = {};
            elseif nargin == 5
                args = {varargin{:}, 'VariableNames', ...
                    {'poly_n', 'poly_coef', 'pwl_n', 'pwl_qty', 'pwl_cost'}};
            else
                error('mp.cost_table: constructor must be called with 0 or 5 arguments.')
            end
            obj@mp_table_subclass(args{:});
        end

        function p = poly_params(obj, idx, pu_base)
            % ::
            %
            %   p = poly_params(obj, idx, pu_base)
            %
            % Inputs:
            %   obj (mp.cost_table) : the cost table
            %   idx: (integer) : index vector of rows of interest, empty
            %       for all rows
            %   pu_base (double) : base used to scale quantities to per unit
            %
            % Outputs:
            %   p (struct) : polynomial cost parameters, struct with fields:
            %
            %     - ``have_quad_cost`` - true if any polynmial costs have
            %       order quadratic or less
            %     - ``i0`` - row indices for constant costs
            %     - ``i1`` - row indices for linear costs
            %     - ``i2`` - row indices for quadratic costs
            %     - ``i3`` - row indices for order 3 or higher costs
            %     - ``k`` - constant term for all quadratic and lower order costs
            %     - ``c`` - linear term for all quadratic and lower order costs
            %     - ``Q`` - quadratic term for all quadratic and lower order costs
            %
            % Implementation in mp.cost_table_utils.poly_params.

            p = mp.cost_table_utils.poly_params(obj, idx, pu_base);
        end

        function p = pwl_params(obj, idx, pu_base, varargin)
            % ::
            %
            %   p = pwl_params(obj, idx, pu_base)
            %   p = pwl_params(obj, idx, pu_base, ng, dc)
            %
            % Inputs:
            %   obj (mp.cost_table) : the cost table
            %   idx: (integer) : index vector of rows of interest, empty
            %       for all rows
            %   pu_base (double) : base used to scale quantities to per unit
            %   ng (integer) : number of units, default is # of rows in cost
            %   dc (boolean) : true if DC formulation (ng variables),
            %       otherwise AC formulation (2*ng variables), default is 1
            %
            % Outputs:
            %   p (struct) : piecewise linear cost parameters, struct with fields:
            %
            %     - ``n`` - number of piecewise linear costs
            %     - ``i`` - row indices for piecewise linear costs
            %     - ``A`` - constraint coefficient matrix for CCV formulation
            %     - ``b`` - constraint RHS vector for CCV formulation
            %
            % Implementation in mp.cost_table_utils.pwl_params.

            p = mp.cost_table_utils.pwl_params(obj, idx, pu_base, varargin{:});
        end

        function maxc = max_pwl_cost(obj)
            % ::
            %
            %   maxc = max_pwl_cost(obj)
            %
            % Input:
            %   obj (mp.cost_table) : the cost table
            %
            % Output:
            %   maxc (double) : maximum cost component of all breakpoints
            %       used to specify piecewise linear costs
            %
            % Implementation in mp.cost_table_utils.max_pwl_cost.

            maxc = mp.cost_table_utils.max_pwl_cost(obj);
        end
    end     %% methods

    methods (Static)
        function [f, df, d2f] = poly_cost_fcn(xx, x_scale, ccm, idx)
            % ::
            %
            %   f = mp.cost_table.poly_cost_fcn(xx, x_scale, ccm, idx)
            %   [f, df] = mp.cost_table.poly_cost_fcn(...)
            %   [f, df, d2f] = mp.cost_table.poly_cost_fcn(...)
            %
            % Evaluates the sum of a set of polynomial cost functions
            % :math:`f(x) = \sum_{i \in I}{f_i(x_i)}`, and optionally the
            % gradient and Hessian.
            %
            % Inputs:
            %   xx (single element cell array of double) : first element is
            %       a vector of the pre-scaled quantities :math:`x/\alpha`
            %       used to compute the costs
            %   x_scale (double) : scalar :math:`\alpha` used to scale the
            %       quantity value before evaluating the polynomial cost
            %   ccm (double) : cost coefficient matrix, element *(i,j)* is
            %       the coefficient of the *(j-1)* order term for cost *i*
            %   idx (integer) : index vector of subset :math:`I` of rows of
            %       ``xx{1}`` and ``ccm`` of interest
            %
            % Outputs:
            %   f (double) : value of cost function :math:`f(x)`
            %   df (vector of double) : (optional) gradient of cost function
            %   d2f (matrix of double) : (optional) Hessian of cost function

            x = xx{1}(idx) * x_scale;
            n = length(xx{1});

            %%----- evaluate cost function -----
            f = sum( mp.cost_table.eval_poly_fcn(ccm(idx, :), x) );

            %%----- evaluate cost gradient -----
            if nargout > 1
                %% coefficients of 1st derivative
                cp = mp.cost_table.diff_poly_fcn(ccm(idx, :));
                df = zeros(n, 1);
                df(idx) = x_scale * mp.cost_table.eval_poly_fcn(cp, x);

                %% ---- evaluate cost Hessian -----
                if nargout > 2
                    %% coefficients of 2nd derivative
                    cpp = mp.cost_table.diff_poly_fcn(cp);
                    d2f = sparse(idx, idx, ...
                        x_scale^2 * mp.cost_table.eval_poly_fcn(cpp, x), n, n);
                end
            end
        end

        function f = eval_poly_fcn(c, x)
            % ::
            %
            %   f = mp.cost_table.eval_poly_fcn(c, x)
            %
            % Evaluate a vector of polynomial functions, where ...
            % ::
            %
            %   f = c(:,1) + c(:,2) .* x + c(:,3) .* x^2 + ...
            %
            % Inputs:
            %   c (matrix of double) : coefficient matrix, element *(i,j)*
            %       is the coefficient of the *(j-1)* order term for *i*-th
            %       element of *f*
            %   x (vector of double) : vector of input values
            %
            % Outputs:
            %   f (vector of double) : value of functions

            if isempty(c)
                f = zeros(size(x));
            else
                f = c(:, 1);        %% constant term
                for k = 2:size(c, 2)
                    f = f + c(:, k) .* x .^ (k-1);
                end
            end
        end

        function c = diff_poly_fcn(c)
            % ::
            %
            %   c = mp.cost_table.diff_poly_fcn(c)
            %
            % Compute the coefficient matrix for the derivatives of a set
            % of polynomial functions from the coefficients of the functions.
            %
            % Inputs:
            %   c (matrix of double) : coefficient matrix for the functions,
            %       element *(i,j)* is the coefficient of the *(j-1)* order
            %       term of the *i*-th function
            %
            % Outputs:
            %   c (matrix of double) : coefficient matrix for the derivatives
            %       of the functions, element *(i,j)* is the coefficient of
            %       the *(j-1)* order term of the derivative of the *i*-th
            %       function

            n = size(c, 2);     %% number of coefficients (cols in c)
            if n >= 2
                c = c(:, 2:n);
            else
                c = zeros(size(c, 1), 1);
            end
            for k = 2:n-1
                c(:, k) = k * c(:, k);
            end
        end
    end     %% methods (Static)
end         %% classdef
