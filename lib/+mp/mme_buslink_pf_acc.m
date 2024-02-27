classdef mme_buslink_pf_acc < mp.mme_buslink_pf_ac
% mp.mme_buslink_pf_acc - Math model element for 1-to-3-phase buslink for AC cartesian voltage PF/CPF.
%
% Math model element class for 1-to-3-phase buslink elements for AC cartesian
% power flow and CPF problems.
%
% Implements methods for adding constraints to match voltages across each
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
            vs_va = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vr_pv', 'Vr3_pv', 'Vr3_pv', 'Vr3_pv', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq', ...
                                    'Vi_pv', 'Vi3_pv', 'Vi3_pv', 'Vi3_pv'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}, ...
                                    {}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            vs_vm = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );

            fcn_va = @(xx)obj.pf_va_fcn(nme, xx, [A_va_pq A_va_pv], b_va);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, [], vs_va);

            fcn_vm = @(xx)obj.pf_vm_fcn(nme, xx, A_vm, b_vm);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, [], vs_vm);

            %% call parent
            add_constraints@mp.mme_buslink_pf_ac(obj, mm, nm, dm, mpopt);
        end

        function [g, dg] = pf_va_fcn(obj, nme, xx, A, b)
            %

            %% unpack data
            vr = vertcat(xx{1:8});
            vi = vertcat(xx{9:16});

            if nargout > 1
                [va, dva] = nme.va_fcn({vr, vi}, [], 0);
                dg = A * dva;
            else
                va = nme.va_fcn({vr, vi}, [], 0);
            end
            g = A * va - b;
        end

        function [g, dg] = pf_vm_fcn(obj, nme, xx, A, b)
            %

            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                n = length(vr);
                [vm, dvm] = nme.vm2_fcn({vr, vi}, [], 0);
                vm_inv = 1./vm;
                dg = 0.5 * A * dvm * spdiags([vm_inv; vm_inv], 0, 2*n, 2*n);
            else
                vm = nme.vm2_fcn({vr, vi}, [], 0);
            end
            g = A * sqrt(vm) - b;
        end
    end     %% methods
end         %% classdef
