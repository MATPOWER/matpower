classdef (Abstract) mm_shared_pfcpf_ac < mp.mm_shared_pfcpf
% mp.mm_shared_pfcpf_ac - Mixin class for AC PF/CPF **math model** objects.
%
% An abstract mixin class inherited by all AC power flow (PF) and continuation
% power flow (CPF) **math model** objects.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function obj = add_system_varset_pf(obj, nm, vvar, typ)
            %

            ad = obj.aux_data;
            st = nm.(vvar);
            d = st.data;
            mmx_i1 = obj.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                ii = ad.node_type_by_elm(k).(typ);
                nii = length(ii);
                if isempty(idx)
                    obj.add_var([name '_' typ], nii, d.v0.(name)(ii), d.vl.(name)(ii), d.vu.(name)(ii));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_' typ], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {ii}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_' typ], idx, nii, v0, vl, vu);
                end
            end
            mmx_iN = obj.var.N;
            if ad.(['n' typ])
                obj.aux_data.var_map{end+1} = ...
                    {vvar, [], [], ad.(typ), mmx_i1, mmx_iN, []};
            end
        end

        function z_ = update_z(obj, nm, v_, z_, ad, Sinj, idx)
            % update_z - Update/allocate active/reactive injections at slack/PV nodes.
            %
            % Update/allocate slack know active power injections and slack/PV
            % node reactive power injections.

            rpv = [ad.ref; ad.pv];      %% slack and PV nodes
            if nargin < 6 || isempty(Sinj)
                %% compute power injection at slack/PV nodes
                idx = find(any(nm.C(rpv, :), 1));   %% ports connected to slack/PV nodes
                Sinj = nm.port_inj_power([v_; z_], 1, idx);
            end
            Sref = nm.C(ad.ref, idx) * Sinj;
            Spv  = nm.C(ad.pv,  idx) * Sinj;
            Qpv = imag(Spv);

            %%-----  active power at slack nodes  -----
            %% coefficient matrix for power injection states
            CC = nm.C * nm.get_params([], 'N') * nm.D';

            %% coefficient matrix for power injection states for slack nodes
            CCref = CC(ad.ref, :);
            jr = find(any(CCref, 1));   %% indices of corresponding states

            %% active power injections at slack nodes
            Pref = real(Sref);

            %% allocate active power at slack nodes to 1st direct inj state
            %% find all z (except first one) with direct injection at each
            %% slack node
            [i, j] = find(CCref);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            %% all included states that are not 1st at their node
            jn = unique(ij(~ismember(1:length(i), k1), 2));

            %% if we have extra states (more than 1) for any node(s)
            if ~isempty(jn)
                %% augment update equation CC * (z - zprev) = -Pref with
                %% additional rows to force these states to remain fixed
                I = speye(nm.nz);
                CCref = [CCref; I(jn, :)];
                Pref = [Pref; zeros(length(jn), 1)];
            end

            %% update z for active injections at slack nodes
            z_(jr) = z_(jr) - CCref(:, jr) \ Pref;

            %%-----  reactive power at slack/PV nodes  -----
            %% coefficient matrix for power injection states for slack/PV nodes
            CCrpv = CC(rpv, :);
            jrpv = find(any(CCrpv, 1));
            Qrpv = [imag(Sref); Qpv] - CCrpv * imag(z_);

            %% find all z with direct injection at each slack/PV node
            [i, j] = find(CCrpv);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            % j1 = ij(k1, 2);     %% indices of states that are 1st at their node
            kn = find(~ismember(1:length(i), k1));  %% indices of entries that are not first state for corresponding node
            %% all included states that are not 1st at their node
            jn = unique(ij(kn, 2));
            in = unique(ij(kn, 1));     %% nodes with multiple states

            %% if we have extra states (more than 1) for some node(s)
            if ~isempty(jn)
                %% find ranges for relevant state vars to allocate reactive
                %% power proportional to ranges
                [~, mn, mx] = nm.params_var('zi');

                %% This code is currently not designed to handle injections
                %% from states that affect more than a single node, so we
                %% throw a warning if we find such a case
                if any(sum(CCrpv ~= 0) > 1)
                    k = find(sum(CCrpv ~= 0) > 1);
                    warning('update_z:multiple_nodes', ...
                        'mp.mm_shared_pfcpf_ac/update_z: unable to distribute reactive power due to z var %d affecting multiple nodes.', k(1));
                end

                %% define a numerical proxy to replace +/- Inf limits
                %% start by setting M equal to average injection at node
                %% CCrpv' * Qrpv is total Qrpv at each corresponding state
                %% CCrpv' * sum(CCrpv, 2) is the number of injections at node
                M = abs((CCrpv' * Qrpv) ./ (CCrpv' * sum(CCrpv, 2)));
                %% add abs value of upper and lower bound to avg nodal injection
                M(~isinf(mx)) = M(~isinf(mx)) + abs(mx(~isinf(mx)));
                M(~isinf(mn)) = M(~isinf(mn)) + abs(mn(~isinf(mn)));
                %% set M for each state to sum over all states at same node
                M = CC' * CC * M;
                %% replace +/- Inf limits with proxy +/- M
                mn(mn ==  Inf) =  M(mn ==  Inf);
                mn(mn == -Inf) = -M(mn == -Inf);
                mx(mx ==  Inf) =  M(mx ==  Inf);
                mx(mx == -Inf) = -M(mx == -Inf);

                %% find indices of states with largest range at its node
                r = mx - mn;    %% size of range for each state
                [rmax, j] = max(abs(CCrpv) * spdiags(r, 0, nm.nz, nm.nz), [], 2);

                %% set ranges to 1 for states at nodes where all ranges are 0
                %% (results in equal limit violations)
                %% find nodes where all corresponding ranges are zero
                i0 = find(abs(rmax) < 10*eps);
                if ~isempty(i0)
                    rmax(i0) = 1;           %% set these ranges to one
                    j0 = find(any(CCrpv(i0, :), 1));    %% corresponding states
                    r(j0) = abs(CCrpv(i0, j0))' * rmax(i0); %% apply to all states at same node
                end

                %% augment update equation ...
                %%  CCrpv * (z - zprev + 1j * imag(zprev)) = -1j * Qrpv
                %% with additional rows to force these states to allocate in
                %% proportion to min-max range, i.e. all states k at node have
                %%  q_k = -qmin_k - r_k * lam, where lam is same for all k at node
                %% we solve for lam using eqn with largest r_k, then
                %% substitute in other equations to get the set of constraints
                %% to add
                R = sparse(nm.nz, nm.nz);
                b = zeros(nm.nz, 1);
                jn = [];    %% initialize list of states that do not have max range at their node
                for i = 1:length(in)    %% for all nodes with multiple states
                    ii = in(i);         %% node ii, with multiple states
                    jj = ij(ij(:, 1) == ii, 2); %% states at node ii
                    rr = r(jj) / r(j(ii));  %% range of each z / max z range
                    R(jj, j(ii)) = rr;  %% rows for states @ node ii, col of state w/max range
                    b(jj) = b(jj) + rr * mn(j(ii)) - mn(jj);
                    jj(jj == j(ii)) = [];
                    jn = [jn; jj];  %% add states at node i that do not have max range
                end
                A = speye(nm.nz) - R;
                Qrpv = [Qrpv; b(jn)];
                CCrpv = [CCrpv; A(jn, :)];
            end

            %% update z for reactive injections at slack/PV nodes
            z0 = z_(jrpv);
            z_(jrpv) = z0 - 1j * (CCrpv(:, jrpv) \ Qrpv + imag(z0));
        end
    end     %% methods
end         %% classdef
