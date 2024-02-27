classdef mme_buslink_pf_acp < mp.mme_buslink_pf_ac
% mp.mme_buslink_pf_acp - Math model element for 1-to-3-phase buslink for AC polar voltage PF/CPF.
%
% Math model element class for 1-to-3-phase buslink elements for AC polar
% power flow and CPF problems.
%
% Implements method for adding constraints to match voltages across each
% buslink.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            nme = obj.network_model_element(nm);

            %% add constraints for matching
            %%  voltage angles at pv and pq nodes
            %%  voltage magnitudes at pq nodes
            [A_va_pq, A_va_pv, b_va, A_vm, b_vm] = obj.voltage_constraints(nme, mm.aux_data);

            %% prep variable set structs
            vs_va = struct('name', {'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv', ...
                                    'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            vs_vm = struct('name', {'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}} );

            %% add constraints
            mm.add_lin_constraint('buslink_va', [A_va_pv A_va_pq], b_va, b_va, vs_va);
            mm.add_lin_constraint('buslink_vm', A_vm, b_vm, b_vm, vs_vm);

            %% call parent
            add_constraints@mp.mme_buslink_pf_ac(obj, mm, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
