classdef cost_table_utils
% mp.cost_table_utils - Static methods for mp.cost_table.
%
% Contains the implementation of some methods that would ideally belong
% in mp.cost_table.
%
% Within classes that inherit from mp_table_subclass, such as mp.cost_table,
% any subscripting to access the elements of the table must be done through
% explicit calls to the table's :meth:`subsref` and :meth:`subsasgn` methods.
% That is, the normal table subscripting syntax will not work, so working
% with the table becomes extremely cumbersome.
%
% This purpose of this class is to provide the implementation for
% mp.cost_table methods that **do** allow access to that table via normal
% table subscripting syntax.
%
% mp.cost_table_util Methods:
%   * poly_params - create struct of polynomial parameters from mp.cost_table
%   * pwl_params - create struct of piecewise linear parameters from mp.cost_table
%   * max_pwl_cost - get maximum cost component used to specify pwl costs
%
% See also mp.cost_table.

%   MATPOWER
%   Copyright (c) 2023-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods (Static)
        function p = poly_params(cost, idx, pu_base)
            % ::
            %
            %   p = mp.cost_table_utils.poly_params(cost, idx, pu_base)
            %
            % Implementation for mp.cost_table.poly_params. See
            % mp.cost_table.poly_params for details.

            if ~isempty(idx)
                cost = cost(idx, :);
            end

            ng = size(cost, 1);
            have_quad_cost = 0;
            kg = []; cg = []; Qg = [];
            i0 = find(cost.poly_n == 1);    %% constant
            i1 = find(cost.poly_n == 2);    %% linear
            i2 = find(cost.poly_n == 3);    %% quadratic
            i3 = find(cost.poly_n > 3);     %% cubic or greater
            if ~isempty(i2) || ~isempty(i1) || ~isempty(i0)
                have_quad_cost = 1;
                kg = zeros(ng, 1);
                cg = zeros(ng, 1);
                if ~isempty(i2)
                    Qg = zeros(ng, 1);
                    Qg(i2) = 2 * cost.poly_coef(i2, 3) * pu_base^2;
                    cg(i2) = cg(i2) + cost.poly_coef(i2, 2) * pu_base;
                    kg(i2) = kg(i2) + cost.poly_coef(i2, 1);
                end
                if ~isempty(i1)
                    cg(i1) = cg(i1) + cost.poly_coef(i1, 2) * pu_base;
                    kg(i1) = kg(i1) + cost.poly_coef(i1, 1);
                end
                if ~isempty(i0)
                    kg(i0) = kg(i0) + cost.poly_coef(i0, 1);
                end
            end
            p = struct( ...
                    'have_quad_cost', have_quad_cost, ...
                    'i0', i0, ...
                    'i1', i1, ...
                    'i2', i2, ...
                    'i3', i3, ...
                    'k', kg, ...
                    'c', cg, ...
                    'Q', Qg ...
                );
        end

        function p = pwl_params(cost, idx, pu_base, ng, dc)
            % ::
            %
            %   p = mp.cost_table_utils.pwl_params(cost, idx, pu_base)
            %   p = mp.cost_table_utils.pwl_params(cost, idx, pu_base, ng, dc)
            %
            % Implementation for mp.cost_table.pwl_params. See
            % mp.cost_table.pwl_params for details.

            if ~isempty(idx)
                cost = cost(idx, :);
            end

            if nargin < 4
                ng = size(cost, 1);
                dc = 1;
            end

            ipwl = find(cost.pwl_n);    %% piece-wise linear costs
            ny = size(ipwl, 1);     %% number of piece-wise linear cost vars
            if dc
                nq = 0;     %% number of qg variables
                q1 = [];
            else
                nq = ng;    %% number of qg variables
                q1 = 1+ng;
            end

            %% from makeAy()
            ybas = ng+nq;
            if ny == 0
                Ay = sparse([], [], [], 0, ybas+ny, 0);
                by = [];
            else
                %% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments,
                %% then the corresponding constraint on pg (or qg) and Y is
                %%                                             c(i+1) - c(i)
                %%  Y   >=   c(i) + m * (pg - p(i)),      m = ---------------
                %%                                             p(i+1) - p(i)
                %%
                %% this becomes   m * pg - Y   <=   m*p(i) - c(i)

                %% form constraint matrix
                m = sum(cost.pwl_n(ipwl));  %% total number of cost points
                Ay = sparse([], [], [], m-ny, ybas+ny, 2*(m-ny));
                by = [];
                k = 1;
                j = 1;
                for i=ipwl'
                    ns = cost.pwl_n(i); %% # of cost points; segments = ns-1
                    p = cost.pwl_qty(i, 1:ns) / pu_base;
                    c = cost.pwl_cost(i, 1:ns);
                    m = diff(c) ./ diff(p);         %% slopes for pg (or qg)
                    if any(diff(p) == 0)
                        fprintf('mp.cost_table_utils/pwl_params: bad qty data in row %i of cost matrix\n', i);
                    end
                    b = m .* p(1:ns-1) - c(1:ns-1); %% and rhs
                    by = [by;  b'];
                    if i > ng
                        sidx = q1 + (i-ng) - 1;     %% for qg cost
                    else
                        sidx = i;                   %% for pg cost
                    end
                    Ay(k:k+ns-2, sidx) = m';
                    Ay(k:k+ns-2, ybas+j) = -ones(ns-1,1);
                    k = k + ns - 1;
                    j = j + 1;
                end
            end

            p = struct('n', ny, 'i', ipwl, 'A', Ay, 'b', by);
        end

        function maxc = max_pwl_cost(cost)
            % ::
            %
            %   maxc = mp.cost_table_utils.max_pwl_cost(cost)
            %
            % Implementation for mp.cost_table.max_pwl_cost. See
            % mp.cost_table.max_pwl_cost for details.

            maxc = max(max(cost.pwl_cost));
        end
    end     %% methods (Static)
end         %% classdef
