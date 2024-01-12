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

        function p = poly_params(cost, idx, pu_base)
            % ::
            %
            %   p = poly_params(cost, idx, pu_base)
            %
            % Inputs:
            %   cost (mp.cost_table) : the cost table
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

            p = mp.cost_table_utils.poly_params(cost, idx, pu_base);
        end

        function p = pwl_params(cost, idx, pu_base, varargin)
            % ::
            %
            %   p = pwl_params(cost, idx, pu_base)
            %   p = pwl_params(cost, idx, pu_base, ng, dc)
            %
            % Inputs:
            %   cost (mp.cost_table) : the cost table
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

            p = mp.cost_table_utils.pwl_params(cost, idx, pu_base, varargin{:});
        end
    end     %% methods
end         %% classdef
