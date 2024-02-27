classdef (Abstract) mme_legacy_dcline_opf < mp.mme_legacy_dcline
% mp.mme_legacy_dcline_opf - Math model element abstract base class for legacy DC line for OPF.
%
% Math model element abstract base class for legacy DC line elements for OPF
% problems.
%
% Implements methods to add costs, including piecewise linear cost variables,
% and to form an interior initial point for cost variables.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % struct for cost parameters with fields:
        %
        %   - ``poly`` - polynomial costs for active power, struct with
        %     fields:
        %
        %       - ``have_quad_cost``
        %       - ``i0``, ``i1``, ``i2``, ``i3``
        %       - ``k``, ``c``, ``Q``
        %   - ``pwl`` - piecewise linear costs for actve power, struct
        %     with fields:
        %
        %       - ``n``, ``i``, ``A``, ``b``
        cost
    end

    methods
        function build_cost_params(obj, dm)
            %
            dme = obj.data_model_element(dm);
            obj.cost = dme.build_cost_params(dm);
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
            %

            %% collect/construct all legacy DC line cost parameters
            obj.build_cost_params(dm);

            %% piecewise linear costs
            if ~isempty(obj.cost) && obj.cost.pwl.n
                mm.add_var('ydc', obj.cost.pwl.n);
            end
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %

            %% add constraint on active flows and ends of DC line
            dme = obj.data_model_element(dm);
            A = [spdiags(dme.loss1 - 1, 0, dme.n, dme.n) -speye(dme.n, dme.n)];
            b = -dme.loss0;
            mm.add_lin_constraint('dcline_legacy', A, b, b, {'Pdcf', 'Pdct'});

            %% piecewise linear costs
            if ~isempty(obj.cost) && obj.cost.pwl.n
                mm.add_lin_constraint('ycondc', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pdcf', 'ydc'});
            end

            %% call parent
            add_constraints@mp.mme_legacy_dcline(obj, mm, nm, dm, mpopt);
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            %

            if ~isempty(obj.cost)
                %% (quadratic) polynomial costs on Pdcf
                if obj.cost.poly.have_quad_cost
                    mm.add_quad_cost('polPdcf', obj.cost.poly.Q, obj.cost.poly.c, obj.cost.poly.k, {'Pdcf'});
                end
    
                %% (order 3 and higher) polynomial costs on Pg
                if ~isempty(obj.cost.poly.i3)
                    dme = obj.data_model_element(dm);
                    cost_Pdcf = @(xx)mp.cost_table.poly_cost_fcn( ...
                        xx, dm.base_mva, dme.tab.cost.poly_coef(dme.on, :), ...
                        obj.cost.poly.i3);
                    mm.add_nln_cost('polPdcf', 1, cost_Pdcf, {'Pdcf'});
                end
    
                %% piecewise linear costs
                if obj.cost.pwl.n
                    mm.add_quad_cost('pwldc', [], ones(obj.cost.pwl.n, 1), 0, {'ydc'});
                end
            end
        end

        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %

            %% set gen cost variables to something feasible
            if ~isempty(obj.cost) && obj.cost.pwl.n > 0
                vv = mm.get_idx();
                dme = obj.data_model_element(dm);
                maxc = max_pwl_cost(dme.tab.cost);
                x0(vv.i1.ydc:vv.iN.ydc) = maxc + 0.1 * abs(maxc);
            end
        end
    end     %% methods
end         %% classdef
