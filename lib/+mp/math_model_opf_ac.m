classdef (Abstract) math_model_opf_ac < mp.math_model_opf
% mp.math_model_opf_ac - Abstract base class for AC OPF **math model** objects.
%
% Provide implementation of nodal current and power balance functions
% and their derivatives, and setup of solver options.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function [g, dg] = nodal_current_balance_fcn(obj, x, nm)
            %

            x_ = obj.convert_x_m2n(x, nm);
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = nm.nodal_complex_current_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% Re{I} mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Im{I} mismatch w.r.t v1, v2, zr, zi
            else
                G = nm.nodal_complex_current_balance(x_);
            end
            g = [ real(G);              %% real current mismatch
                  imag(G) ];            %% imaginary current mismatch
        end

        function [g, dg] = nodal_power_balance_fcn(obj, x, nm)
            %

            x_ = obj.convert_x_m2n(x, nm);
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = nm.nodal_complex_power_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = nm.nodal_complex_power_balance(x_);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end

        function d2G = nodal_current_balance_hess(obj, x, lam, nm)
            %

            x_ = obj.convert_x_m2n(x, nm);
            nlam = length(lam) / 2;
            lamIr = lam(1:nlam);
            lamIi = lam((1:nlam)+nlam);

            d2Gr = nm.nodal_complex_current_balance_hess(x_, lamIr);
            d2Gi = nm.nodal_complex_current_balance_hess(x_, lamIi);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function d2G = nodal_power_balance_hess(obj, x, lam, nm)
            %

            x_ = obj.convert_x_m2n(x, nm);
            nlam = length(lam) / 2;
            lam_p = lam(1:nlam);
            lam_q = lam((1:nlam)+nlam);

            d2Gr = nm.nodal_complex_power_balance_hess(x_, lam_p);
            d2Gi = nm.nodal_complex_power_balance_hess(x_, lam_q);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %

            opt = mpopt2nlpopt(mpopt, obj.problem_type());

            if mpopt.opf.start < 2      %% initialize interior point
                opt.x0 = obj.interior_x0(obj, nm, dm);
            end
        end
    end     %% methods
end         %% classdef
