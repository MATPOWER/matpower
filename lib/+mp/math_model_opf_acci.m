classdef math_model_opf_acci < mp.math_model_opf_acc
% mp.math_model_opf_acci - OPF **math model** for AC-cartesian-current formulation.
%
% Implements formulation-specific and OPF-specific node balance constraint
% and node balance price methods.

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
        function tag = form_tag(obj)
            %

            tag = 'acci';
        end

        function name = form_name(obj)
            %

            name = 'AC-cartesian-current';
        end

        function add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            %% power balance constraints
            nn = nm.node.N;             %% number of nodes
            fcn_mis = @(x)obj.nodal_current_balance_fcn(x, nm);
            hess_mis = @(x, lam)obj.nodal_current_balance_hess(x, lam, nm);
            obj.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lam_p, lam_q] = node_power_balance_prices(obj, nm)
            %

            %% shadow prices on node power balance
            nne = obj.get_idx('nle');

            %% convert current balance shadow prices to equivalent lam_p and lam_q
            %% P + jQ = (Vr + jVi) * (M - jN)
            %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
            %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
            %% lam_p = df/dP = df/dM * dM/dP + df/dN + dN/dP
            %% lam_q = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
            V = nm.soln.v;
            lambda = obj.soln.lambda;
            VV = V ./ (V .* conj(V));   %% V / vm^2
            VVr = real(VV);
            VVi = imag(VV);
            lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
            lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
            lam_p = (VVr.*lamM + VVi.*lamN);
            lam_q = (VVi.*lamM - VVr.*lamN);
        end
    end     %% methods
end         %% classdef
