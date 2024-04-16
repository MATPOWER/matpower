classdef (Abstract) mm_shared_pfcpf_dc < mp.mm_shared_pfcpf
% mp.mm_shared_pfcpf_dc - Mixin class for DC power flow (PF) **math model** objects.
%
% An abstract mixin class inherited by DC power flow (PF) **math model**
% objects.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function ad = build_aux_data(obj, nm, dm, mpopt)
            %

            %% call parent
            ad = build_aux_data@mp.mm_shared_pfcpf(obj, nm, dm, mpopt);

            %% get parameters
            [B, K, p] = nm.get_params();

            %% add DC model parameters
            ad.B = nm.C * B * nm.C';
            ad.Pbus = -(nm.C * K * nm.D' * ad.z + nm.C * p);
        end

        function obj = add_system_vars_pf(obj, nm, dm, mpopt)
            %

            %% get model variables
            vvars = nm.model_vvars();

            %% index vectors
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = nm.(vvars{1});
            d = st.data;
            mmx_i1 = obj.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    obj.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp.mm_shared_pfcpf_dc.add_system_vars_pf: handling of indexed sets not implmented here (yet)');
                end
            end
            mmx_iN = obj.var.N;
            if ad.npv || ad.npq
                obj.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], pvq, mmx_i1, mmx_iN, []};
            end
        end

        function [vx, z, x] = convert_x_m2n(obj, mmx, nm, only_v)
            % convert_x_m2n - Convert math model state to network model state.
            % ::
            %
            %   x = mm.pf_convert(mmx, nm)
            %   [v, z] = mm.pf_convert(mmx, nm)
            %   [v, z, x] = mm.pf_convert(mmx, nm)
            %   ... = mm.pf_convert(mmx, nm, only_v)

            %% update v_, z_ from mmx
            nm_vars = obj.update_nm_vars(mmx, nm);
            vx = nm_vars.va;
            z  = nm_vars.z;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z = obj.update_z(nm, vx, z, obj.aux_data);
            end

            %% prepare return values
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function z = update_z(obj, nm, v, z, ad)
            % update_z - Update/allocate slack node active power injections.

            %% coefficient matrix for power injection states at slack bus
            CC = nm.C(ad.ref, :) * nm.get_params([], 'K') * nm.D';
            jr = find(any(CC, 1));  %% indices of corresponding states

            %% power injections at slack nodes
            idx = find(any(nm.C(ad.ref, :), 1));    %% ports connected to slack nodes
            Pref = nm.C(ad.ref, idx) * nm.port_inj_power([v; z], 1, idx);

            %% allocate active power at slack nodes to 1st direct inj state
            %% find all z (except first one) with direct injection at each
            %% slack node
            [i, j] = find(CC);
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
                CC = [CC; I(jn, :)];
                Pref = [Pref; zeros(length(jn), 1)];
            end

            %% update z for active injections at slack nodes
            z(jr) = z(jr) - CC(:, jr) \ Pref;
        end
    end     %% methods
end         %% classdef
