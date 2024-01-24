classdef mme_bus_nld_opf_acps_node_test < mp.mme_bus_opf_acp

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function name = name(obj)
            name = 'bus_nld';
        end

        function x0 = interior_x0(obj, mm, nm, dm, x0)
            varef1 = mm.interior_va(nm, dm);
            vm = obj.interior_vm(mm, nm, dm);
            vv = mm.get_idx();
            x0(vv.i1.(['va_' obj.name]):vv.iN.(['va_' obj.name])) = varef1; %% angles set to 1st ref angle
            x0(vv.i1.(['vm_' obj.name]):vv.iN.(['vm_' obj.name])) = vm;     %% voltage magnitudes
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.(obj.name):nn.iN.(obj.name));

            %% shadow prices on voltage magnitudes
            vv = mm.get_idx('var');
            lambda = mm.soln.lambda;
            vVm = ['vm_' obj.name];
            mu_vm_lb = lambda.lower(vv.i1.(vVm):vv.iN.(vVm));
            mu_vm_ub = lambda.upper(vv.i1.(vVm):vv.iN.(vVm));

            %% shadow prices on node power balance
            [lam_p, lam_q] = mm.node_power_balance_prices(nm);
            lam_p = lam_p(nn.i1.(obj.name):nn.iN.(obj.name));   %% for (obj.name) nodes only
            lam_q = lam_q(nn.i1.(obj.name):nn.iN.(obj.name));   %% for (obj.name) nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = angle(V) * 180/pi;
            dme.tab.vm(dme.on) = abs(V);
            dme.tab.lam_p(dme.on) = lam_p / dm.base_mva;
            dme.tab.lam_q(dme.on) = lam_q / dm.base_mva;
            dme.tab.mu_vm_lb(dme.on) = mu_vm_lb;
            dme.tab.mu_vm_ub(dme.on) = mu_vm_ub;
        end
    end     %% methods
end         %% classdef
