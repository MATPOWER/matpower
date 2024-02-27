classdef nme_bus_dc < mp.nme_bus & mp.form_dc
% mp.nme_bus_dc - Network model element for bus for DC formulations.
%
% Adds voltage variable ``Va`` to the network model and inherits from
% mp.form_dc.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_vvars(obj, nm, dm, idx)
            %
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare angle bounds for ref buses
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);
            k = find(dme.type == mp.NODE_TYPE.REF);
            va_lb(k) = dme.va_start(k);
            va_ub(k) = dme.va_start(k);

            nm.add_var('va', 'Va', nb, dme.va_start, va_lb, va_ub);
        end
    end     %% methods
end         %% classdef
