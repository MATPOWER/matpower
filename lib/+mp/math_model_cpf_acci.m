classdef math_model_cpf_acci < mp.math_model_cpf_acc & mp.mm_shared_pfcpf_acci
% mp.math_model_cpf_acci - CPF **math model** for AC-cartesian-current formulation.
%
% Implements formulation-specific and CPF-specific node balance constraint.

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

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)node_balance_equations_cpf(obj, x, nm);
            obj.add_nln_constraint({'Irmis', 'Iimis', 'Vmis'}, [npvq;npvq;ad.npv], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
