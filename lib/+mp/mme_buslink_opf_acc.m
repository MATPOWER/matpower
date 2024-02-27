classdef mme_buslink_opf_acc < mp.mme_buslink_opf
% mp.mme_buslink_opf_acc - Math model element for 1-to-3-phase buslink for AC cartesian voltage OPF.
%
% Math model element class for 1-to-3-phase buslink elements for AC cartesian
% OPF problems.
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

            %% voltage equality constraints
            [A, b_va, b_vm] = nme.voltage_constraints();
            vs = struct('name', {'Vr', 'Vr3', 'Vr3', 'Vr3', ...
                                 'Vi', 'Vi3', 'Vi3', 'Vi3'}, ...
                        'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}});

            fcn_va = @(xx)obj.va_fcn(nme, xx, A, b_va);
            hess_va = @(xx, lam)obj.va_hess(nme, xx, lam, A);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, hess_va, vs);

            fcn_vm = @(xx)obj.vm2_fcn(nme, xx, A, b_vm);
            hess_vm = @(xx, lam)obj.vm2_hess(nme, xx, lam, A);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, hess_vm, vs);
        end

        function [g, dg] = va_fcn(obj, nme, xx, A, b)
            %

            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                [va, dva] = nme.va_fcn({vr, vi}, [], 0);
                dg = A * dva;
            else
                va = nme.va_fcn({vr, vi}, [], 0);
            end
            g = A * va - b;
        end

        function d2G = va_hess(obj, nme, xx, lam, A)
            %

            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            d2G = nme.va_hess({vr, vi}, A' * lam, []);
        end

        function [g, dg] = vm2_fcn(obj, nme, xx, A, b)
            %

            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                [vm, dvm] = nme.vm2_fcn({vr, vi}, [], 0);
                dg = A * dvm;
            else
                vm = nme.vm2_fcn({vr, vi}, [], 0);
            end
            g = A * vm - b;
        end

        function d2G = vm2_hess(obj, nme, xx, lam, A)
            %

            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            d2G = nme.vm2_hess({vr, vi}, A' * lam, []);
        end
    end     %% methods
end         %% classdef
