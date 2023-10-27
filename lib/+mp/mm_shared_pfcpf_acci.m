classdef (Abstract) mm_shared_pfcpf_acci < mp.mm_shared_pfcpf_acc & mp.mm_shared_pfcpf_ac_i
% mp.mm_shared_pfcpf_acci - Mixin class for AC-cartesian-current PF/CPF **math model** objects.
%
% An abstract mixin class inherited by AC power flow (PF) and continuation
% power flow (CPF) **math model** objects that use a cartesian voltage
% and current balance formuation.

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
            ad = build_aux_data@mp.mm_shared_pfcpf_acc(obj, nm, dm, mpopt);

            %% add data needed for current formulations
            ad = obj.build_aux_data_i(nm, ad);
        end

        function obj = add_system_vars_pf(obj, nm, dm, mpopt)
            %

            %% get model variables
            vvars = nm.model_vvars();

            %% reactive injections
            ad = obj.aux_data;
            v_ = ad.vr + 1j * ad.vi;
            z_ = ad.zr + 1j * ad.zi;
            Qpv = nm.C(ad.pv, :) * imag( nm.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv + ad.zi(ad.zi_idx);
            mmx_i1 = obj.var.N + 1;
            obj.add_var('Qg_pv', ad.npv, Qg_pv);
            mmx_iN = obj.var.N;
            if ad.npv
                obj.aux_data.var_map{end+1} = ...
                    {'zi', [], [], ad.zi_idx, mmx_i1, mmx_iN, []};
            end

            %% voltage real part
            obj.add_system_varset_pf(nm, vvars{1}, 'pq');
            obj.add_system_varset_pf(nm, vvars{1}, 'pv');

            %% voltage imaginary part
            obj.add_system_varset_pf(nm, vvars{2}, 'pq');
            obj.add_system_varset_pf(nm, vvars{2}, 'pv');
        end

        function [f, J] = node_balance_equations(obj, x, nm)
            %

            %% index vectors
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.convert_x_m2n(x, nm, 1);

            %% incidence matrix
            C = nm.C;

            %% Jacobian
            if nargout > 1
                %% get port current injections with derivatives
                [I, dI.vr, dI.vi, dI.zr, dI.zi] = nm.port_inj_current([v_; z_], 1);
                dI.vr = C * dI.vr;
                dI.vi = C * dI.vi;
                dI.zr = C * dI.zr;
                dI.zi = C * dI.zi;
                %% derivatives of voltage magnitudes (for PV buses)
                nn = nm.node.N;
                dV2.vr = sparse(ad.pv, ad.pv, 2*real(v_(ad.pv)), nn, nn);
                dV2.vi = sparse(ad.pv, ad.pv, 2*imag(v_(ad.pv)), nn, nn);
                dV2.zr = sparse(nn, nn);
                dV2.zi = dV2.zr;
                JJ = cell(3, length(ad.var_map));

                for k = 1:length(ad.var_map)
                    m = ad.var_map{k};
                    name = m{1};
                    if ~isempty(m{2})       %% i1:iN
                        i1 = m{2};
                        iN = m{3};
                        JJ{1, k} = real(dI.(name)(pvq, i1:iN));
                        JJ{2, k} = imag(dI.(name)(pvq, i1:iN));
                        JJ{3, k} = dV2.(name)(ad.pv, i1:iN);
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dI.(name)(pvq, :));
                        JJ{2, k} = imag(dI.(name)(pvq, :));
                        JJ{3, k} = dV2.(name)(ad.pv, :);
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dI.(name)(pvq, idx));
                        JJ{2, k} = imag(dI.(name)(pvq, idx));
                        JJ{3, k} = dV2.(name)(ad.pv, idx);
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :}), ...
                             horzcat(JJ{3, :})  );
            else
                %% get port current injections (w/o derivatives)
                I = nm.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.vr(ad.pv).^2 - ad.vi(ad.pv).^2;
            f = [real(II(pvq)); imag(II(pvq)); vmm];
        end
    end     %% methods
end         %% classdef
