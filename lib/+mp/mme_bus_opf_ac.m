classdef (Abstract) mme_bus_opf_ac < mp.mme_bus

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end

    methods
        function vm = interior_vm(obj, mm, nm, dm)
            %% return vm equal to avg of clipped limits
            dme = obj.data_model_element(dm);
            vm_ub = min(dme.vm_ub, 1.5);
            vm_lb = max(dme.vm_lb, 0.5);
            vm = (vm_ub + vm_lb) / 2;
        end
    end     %% methods
end         %% classdef
