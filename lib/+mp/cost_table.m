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
%   the object-oriented ``T.<method>(...)``, to call mp.cost_table methods.
%
% mp.cost_table Methods:
%   * cost_table - construct object
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
% See also mp_table, mp_table_subclass.

%   MATPOWER
%   Copyright (c) 2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = cost_table(poly_n, poly_coef, pwl_n, pwl_qty, pwl_cost)
            % ::
            %
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

            obj.tab = mp_table_subclass(poly_n, poly_coef, pwl_n, pwl_qty, pwl_cost, ...
                'VariableNames', ...
                {'poly_n', 'poly_coef', 'pwl_n', 'pwl_qty', 'pwl_cost'});
        end
    end     %% methods
end         %% classdef
