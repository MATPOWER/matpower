classdef (Abstract) mme_gen_opf < mp.mme_gen
% mp.mme_gen_opf - Math model element abstract base class for generator for OPF.
%
% Math model element abstract base class for generator elements for OPF
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
        %   - ``poly_p`` - polynomial costs for active power, struct
        %     returned by mp.cost_table.poly_params, with fields:
        %
        %       - ``have_quad_cost``
        %       - ``i0``, ``i1``, ``i2``, ``i3``
        %       - ``k``, ``c``, ``Q``
        %   - ``poly_q`` - polynomial costs for reactive power
        %     *(same struct as* ``poly_p`` *)*
        %   - ``pwl`` - piecewise linear costs for actve & reactive struct
        %     returned by mp.cost_table.pwl_params, with fields:
        %
        %       - ``n``, ``i``, ``A``, ``b``
        cost
    end

    methods
        function obj = add_vars(obj, mm, nm, dm, mpopt)
            %

            %% collect/construct all generator cost parameters
            obj.build_cost_params(dm);

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_var('y', obj.cost.pwl.n);
            end
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            %

            %% (quadratic) polynomial costs on Pg
            if obj.cost.poly_p.have_quad_cost
                mm.add_quad_cost('polPg', obj.cost.poly_p.Q, obj.cost.poly_p.c, obj.cost.poly_p.k, {'Pg'});
            end

            %% (order 3 and higher) polynomial costs on Pg
            if ~isempty(obj.cost.poly_p.i3)
                dme = obj.data_model_element(dm);
                cost_Pg = @(xx)mp.cost_table.poly_cost_fcn( ...
                    xx, dm.base_mva, ...
                    dme.tab.cost_pg.poly_coef(dme.on, :), ...
                    obj.cost.poly_p.i3);
                mm.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
            end

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_quad_cost('pwl', [], ones(obj.cost.pwl.n, 1), 0, {'y'});
            end
        end

        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %

            %% set gen cost variables to something feasible
            if obj.cost.pwl.n > 0
                vv = mm.get_idx();
                dme = obj.data_model_element(dm);
                maxgc = dme.max_pwl_gencost();
                x0(vv.i1.y:vv.iN.y) = maxgc + 0.1 * abs(maxgc);
            end
        end
    end     %% methods
end         %% classdef
