classdef (Abstract) math_model_opf_acp < mp.math_model_opf_ac
% mp.math_model_opf_acp - Abstract base class for AC polar OPF **math model** objects.
%
% Provides formulation-specific and OPF-specific subclasses for elements.
%
% Implements convert_x_m2n() to convert from math model state to network
% model state.

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
        %% constructor
        function obj = math_model_opf_acp()
            %

            obj@mp.math_model_opf_ac();
            obj.element_classes = { @mp.mme_bus_opf_acp, @mp.mme_gen_opf_ac, ...
                @mp.mme_load_pf_ac, @mp.mme_branch_opf_acp, @mp.mme_shunt_pf_ac };
        end

        function [vx_, z_, x_] = convert_x_m2n(obj, mmx, nm)
            %

            nm_vars = obj.update_nm_vars(mmx, nm);

            %% convert (real) math model x to (complex) network model x_
            vx_ = nm_vars.vm .* exp(1j * nm_vars.va);
            z_  = nm_vars.zr + 1j * nm_vars.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end
    end     %% methods
end         %% classdef
