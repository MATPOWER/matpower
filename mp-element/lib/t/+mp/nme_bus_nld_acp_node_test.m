classdef nme_bus_nld_acp_node_test < mp.nme_bus_acp

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'bus_nld';
        end

        function obj = add_vvars(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare angle bounds for ref buses
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);
            k = find(dme.type == mp.NODE_TYPE.REF);
            va_lb(k) = dme.va_start(k);
            va_ub(k) = dme.va_start(k);

            nm.add_var('va', ['va_' obj.name], nb, dme.va_start, va_lb, va_ub);
            nm.add_var('vm', ['vm_' obj.name], nb, dme.vm_start, dme.vm_lb, dme.vm_ub);
        end
    end     %% methods
end         %% classdef
