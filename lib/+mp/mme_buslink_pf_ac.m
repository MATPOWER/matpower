classdef (Abstract) mme_buslink_pf_ac < mp.mme_buslink
% mp.mme_buslink_pf_ac - Math model element abstract base class for 1-to-3-phase buslink for AC PF/CPF.
%
% Abstract math model element base class for 1-to-3-phase buslink elements
% for AC power flow and CPF problems.
%
% Implements methods for adding per-phase active and reactive power variables
% and for forming and adding voltage and reactive power constraints.

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
        function obj = add_vars(obj, mm, nm, dm, mpopt)
            %
            nme = obj.network_model_element(nm);

            mm.init_indexed_name('var', 'Plink', {3});
            mm.init_indexed_name('var', 'Qlink', {3});
            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Plink', {p}, nme.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zr', nm.zr.idx.i1.Plink(1), nm.zr.idx.iN.Plink(end), [], mmx_i1, mmx_iN, []};

            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Qlink', {p}, nme.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zi', nm.zi.idx.i1.Qlink(1), nm.zi.idx.iN.Qlink(end), [], mmx_i1, mmx_iN, []};
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            nme = obj.network_model_element(nm);

            %% add Qlink constraints for PV node on 3-phase side

            %% indices of buslinks connected to PV nodes via port 2 (3&4)
            ad = mm.aux_data;
            [~, k, ~] = find(nme.C(ad.pv, nme.nk+1:2*nme.nk));
            n = length(k);
            if n
                I = sparse(1:n, k, 1, n, nme.nk);
                zz = sparse(n, nme.nk);
                A = [I -I zz; I zz -I];
                b = zeros(2*n, 1);
                vs = struct('name', {'Qlink', 'Qlink', 'Qlink'}, ...
                            'idx', {{1}, {2}, {3}} );
                mm.add_lin_constraint('buslink_qlink', A, b, b, vs);
            end
        end

        function [A_va_pq, A_va_pv, b_va, A_vm, b_vm] = voltage_constraints(obj, nme, ad)
            %

            %% form constraint matrices for matching
            %%  voltage angles for pv and pq nodes
            %%  voltage magnitudes for pq nodes
            %% columns of A_va_pv correspond to ...
            %%  'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv',
            %% columns of A_va_pq correspond to ...
            %%  'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'
            %% columns of A_vm correspond to ...
            %%  'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'
            nk = nme.nk;

            %% basic constraint matrix for voltage equality (all nodes)
            [A, b_va, b_vm, Istack] = nme.voltage_constraints();

            %% sub-matrices by node type
            A_ref = A(:, ad.ref);
            A_pv = A(:, ad.pv);
            A_pq = A(:, ad.pq);

            %% voltage angle constraints (all buslinks)
            A_va_pq = A_pq;
            A_va_pv = A_pv;

            %% voltage magnitude constraints (all buslinks)
            A_vm = A_pq;

            %% indices of buslinks connected to REF nodes (fixed va)
            [~, k_va1, ~] = find(nme.C(ad.ref, 1:nk));      %% via port 1
            [~, k_va2, ~] = find(nme.C(ad.ref, nk+1:2*nk)); %% via port 2 (3&4)

            %% adjust RHS of va constraints involving fixed voltage angles
            if ~isempty(k_va1) || ~isempty(k_va2)
                [va, vm] = nme.aux_data_va_vm(ad);
                va1 = nme.C(:, 1:nk)' * va;         %% port 1 voltage angles
                va2 = nme.C(:, nk+1:2*nk)' * va;    %% port 2 voltage angles
                va_adj = zeros(nk, 1);
                va_adj(k_va1) = -va1(k_va1);
                va_adj(k_va2) =  va2(k_va2);
                b_va = b_va + Istack * va_adj;

                %% delete redundate (and currently incorrect) constraints
                %% corresponding to k_va2 and ports 3 and 4
                if ~isempty(k_va2)
                    d_va = [nk+k_va2; 2*nk+k_va2];
                    A_va_pq(d_va, :) = [];
                    A_va_pv(d_va, :) = [];
                    b_va(d_va) = [];
                end
            end

            %% indices of buslinks connected to REF/PV nodes (fixed vm)
            rpv = [ad.ref;ad.pv];   %% indices of REF/PV nodes
            [~, k_vm1, ~] = find(nme.C(rpv, 1:nk));         %% via port 1
            [~, k_vm2, ~] = find(nme.C(rpv, nk+1:2*nk));    %% via port 2 (3&4)

            %% adjust RHS of vm constraints involving fixed voltage magnitudes
            if ~isempty(k_vm1) || ~isempty(k_vm2)
                [va, vm] = nme.aux_data_va_vm(ad);
                vm1 = nme.C(:, 1:nk)' * vm;         %% port 1 voltage magnitudes
                vm2 = nme.C(:, nk+1:2*nk)' * vm;    %% port 2 voltage magnitudes
                vm_adj = zeros(nk, 1);
                vm_adj(k_vm1) = -vm1(k_vm1);
                vm_adj(k_vm2) =  vm2(k_vm2);
                b_vm = b_vm + Istack * vm_adj;

                %% delete redundate (and currently incorrect) constraints
                %% corresponding to k_vm2 and ports 3 and 4
                if ~isempty(k_vm2)
                    d_vm = [nk+k_vm2; 2*nk+k_vm2];
                    A_vm(d_vm, :) = [];
                    b_vm(d_vm) = [];
                end
            end
        end
    end     %% methods
end         %% classdef
