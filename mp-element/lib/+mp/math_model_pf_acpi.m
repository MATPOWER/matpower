classdef math_model_pf_acpi < mp.math_model_pf & mp.mm_shared_pfcpf_acpi
%MP.MATH_MODEL_PF_ACPI  MATPOWER mathematical model for AC power flow (PF) problem.
%   ?
%
%   MP.MATH_MODEL_PF_ACPI ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = math_model_pf_acpi()
            obj@mp.math_model_pf();
            obj.element_classes = { @mp.mme_bus_pf_acp, @mp.mme_gen_pf_ac, ...
                @mp.mme_load_pf_ac, @mp.mme_branch_pf_ac, @mp.mme_shunt_pf_ac };
        end

        function tag = form_tag(obj)
            tag = 'acpi';
        end

        function name = form_name(obj)
            name = 'AC-polar-current';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)node_balance_equations(obj, x, nm);
            obj.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
