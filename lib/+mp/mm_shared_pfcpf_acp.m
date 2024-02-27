classdef (Abstract) mm_shared_pfcpf_acp < mp.mm_shared_pfcpf_ac
% mp.mm_shared_pfcpf_acp - Mixin class for AC polar PF/CPF **math model** objects.
%
% An abstract mixin class inherited by all AC power flow (PF) and continuation
% power flow (CPF) **math model** objects that use a polar voltage
% formuation.

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
        function [vx_, z_, x_] = convert_x_m2n(obj, mmx, nm, only_v)
            % convert_x_m2n - Convert math model state to network model state.
            % ::
            %
            %   x = mm.pf_convert(mmx, nm)
            %   [v, z] = mm.pf_convert(mmx, nm)
            %   [v, z, x] = mm.pf_convert(mmx, nm)
            %   ... = mm.pf_convert(mmx, nm, only_v)

            %% update v_, z_ from mmx
            nm_vars = obj.update_nm_vars(mmx, nm);
            vx_ = nm_vars.vm .* exp(1j * nm_vars.va);
            z_ = nm_vars.zr + 1j * nm_vars.zi;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z_ = obj.update_z(nm, vx_, z_, obj.aux_data);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end
    end     %% methods
end         %% classdef
