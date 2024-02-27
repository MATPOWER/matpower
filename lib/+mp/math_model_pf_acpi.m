classdef math_model_pf_acpi < mp.math_model_pf_ac & mp.mm_shared_pfcpf_acpi
% mp.math_model_pf_acpi - Power flow (PF) **math model** for AC-polar-current formulation.
%
% Implements formulation-specific node balance constraints and inherits
% from formulation-specific class for shared PF/CPF code.

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

            tag = 'acpi';
        end

        function name = form_name(obj)
            %

            name = 'AC-polar-current';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)node_balance_equations(obj, x, nm);
            obj.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
