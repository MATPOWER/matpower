classdef (Abstract) mm_shared_pfcpf_acps < mp.mm_shared_pfcpf_acp
% mp.mm_shared_pfcpf_acps - Mixin class for AC-polar-power PF/CPF **math model** objects.
%
% An abstract mixin class inherited by AC power flow (PF) and continuation
% power flow (CPF) **math model** objects that use a polar voltage
% and power balance formuation.

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
            ad = build_aux_data@mp.mm_shared_pfcpf_acp(obj, nm, dm, mpopt);

            switch mpopt.pf.alg
                case 'GS'
                    ad.Y = nm.C * nm.get_params([], 'Y') * nm.C';
                case 'ZG'
                    pvq = [ad.pv; ad.pq];
                    Y = nm.C * nm.get_params([], 'Y') * nm.C';
                    Y21 = Y(pvq, ad.ref);
                    Y22 = Y(pvq, pvq);
                    [L, U, p, q] = lu(Y22, 'vector');
                    [junk, iq] = sort(q);
                    [ad.Y, ad.Y21, ad.L, ad.U, ad.p, ad.iq] = ...
                        deal(Y, Y21, L, U, p, iq);
            end
        end

        function obj = add_system_vars_pf(obj, nm, dm, mpopt)
            %

            %% get model variables
            vvars = nm.model_vvars();

            %% voltage angles
            obj.add_system_varset_pf(nm, vvars{1}, 'pv');
            obj.add_system_varset_pf(nm, vvars{1}, 'pq');

            %% voltage magnitudes
            obj.add_system_varset_pf(nm, vvars{2}, 'pq');
        end

        function [f, J] = node_balance_equations(obj, x, nm, fdpf)
            %

            %% index vector
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.convert_x_m2n(x, nm, 1);

            %% incidence matrix
            C = nm.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                var_names = cellfun(@(x)x{1}, ad.var_map, 'UniformOutput', false);
                dz = any(strcmp(var_names, 'zr')) || ...
                     any(strcmp(var_names, 'zi'));
                if dz
                    [S, dS.va, dS.vm, dS.zr, dS.zi] = nm.port_inj_power([v_; z_], 1);
                else
                    [S, dS.va, dS.vm] = nm.port_inj_power([v_; z_], 1);
                end
                dS.va = C * dS.va;
                dS.vm = C * dS.vm;
                if dz
                    dS.zr = C * dS.zr;
                    dS.zi = C * dS.zi;
                end
                JJ = cell(2, length(ad.var_map));

                for k = 1:length(ad.var_map)
                    m = ad.var_map{k};
                    name = m{1};
                    if ~isempty(m{2})       %% i1:iN
                        i1 = m{2};
                        iN = m{3};
                        JJ{1, k} = real(dS.(name)(pvq,   i1:iN));
                        JJ{2, k} = imag(dS.(name)(ad.pq, i1:iN));
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dS.(name)(pvq,   :));
                        JJ{2, k} = imag(dS.(name)(ad.pq, :));
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dS.(name)(pvq,   idx));
                        JJ{2, k} = imag(dS.(name)(ad.pq, idx));
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :})  );
            else
                %% get port power injections (w/o derivatives)
                S = nm.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance
            if nargin > 3 && fdpf
                SS = C * S ./ abs(v_);  %% for fast-decoupled formulation
            else
                SS = C * S;
            end
            f = [real(SS(pvq)); imag(SS(ad.pq))];
        end
    end     %% methods
end         %% classdef
