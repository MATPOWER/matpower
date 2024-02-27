classdef (Abstract) nme_buslink < mp.nm_element %& mp.form_ac
% mp.nme_buslink - Network model element abstract base class for 1-to-3-phase buslink.
%
% Implements the network model element for 1-to-3-phase buslink elements,
% with 4 ports and 3 non-voltage states per buslink.
%
% Adds non-voltage state variables ``Plink`` and ``Qlink`` to the network
% model, builds the parameter :math:`\NN`, and constructs voltage constraints.

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
        function name = name(obj)
            %
            name = 'buslink';
        end

        function np = np(obj)
            %
            np = 4;     %% this is a 4 port element
        end

        function nz = nz(obj)
            %
            nz = 3;     %% 3 complex non-voltage states per element
        end

        function obj = add_zvars(obj, nm, dm, idx)
            %
            p = idx{1};
            ng = obj.nk;

            if p == 1
                nm.init_indexed_name('zr', 'Plink', {obj.nz});
                nm.init_indexed_name('zi', 'Qlink', {obj.nz});
            end
            nm.add_var('zr', 'Plink', {p}, obj.nk, 0, -Inf, Inf);
            nm.add_var('zi', 'Qlink', {p}, obj.nk, 0, -Inf, Inf);
        end

        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nm_element(obj, nm, dm);    %% call parent
            I = (dm.base_kva / dm.base_mva / 1000) * speye(obj.nk);
            obj.N = [ repmat(I, 1, obj.nz);
                     -speye(obj.nk * obj.nz) ];
        end

        function [A, b_va, b_vm, Istack_] = voltage_constraints(obj)
            %

            %% form constraint matrices for matching voltages
            nk = obj.nk;

            %% basic constraint matrix for voltage equality (all nodes)
            Istack = repmat(speye(nk), obj.nz, 1);  %% stacked identities
            A = [ Istack -speye(nk * obj.nz) ] * obj.C';
            ang120 = 2*pi/3*ones(nk, 1);

            %% RHS
            b_va = [zeros(nk, 1); ang120; -ang120]; %% angle constraints
            b_vm = zeros(nk*obj.nz, 1);             %% magnitude constraints

            if nargout > 3
                Istack_ = Istack;
            end
        end
    end     %% methods
end         %% classdef
