classdef nme_bus3p_acp < mp.nme_bus3p & mp.form_acp
% mp.nme_bus3p_acp - Network model element for 3-phase bus, AC polar voltage formulation.
%
% Adds voltage variables ``Va3`` and ``Vm3`` to the network model and inherits
% from mp.form_acp.

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
        function obj = add_vvars(obj, nm, dm, idx)
            %
            dme = obj.data_model_element(dm);
            nb = obj.nk;
            p = idx{1};

            %% prepare angle bounds for ref buses
            ref = dme.type == mp.NODE_TYPE.REF;
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);
            vm_start = dme.(sprintf('vm%d_start', p));
            va_start = dme.(sprintf('va%d_start', p));
            va1_lb(ref) = va_start(ref);
            va1_ub(ref) = va_start(ref);

            if p == 1
                nm.init_indexed_name('va', 'Va3', {obj.nn});
                nm.init_indexed_name('vm', 'Vm3', {obj.nn});
            end
            nm.add_var('va', 'Va3', {p}, nb, va_start, va_lb, va_ub);
            nm.add_var('vm', 'Vm3', {p}, nb, vm_start, 0, 2);
        end
    end     %% methods
end         %% classdef
