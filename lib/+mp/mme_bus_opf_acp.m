classdef mme_bus_opf_acp < mp.mme_bus_opf_ac
% mp.mme_bus_opf_acp - Math model element for bus for AC polar voltage OPF.
%
% Math model element class for bus elements for AC polar voltage OPF
% problems.
%
% Implements methods for forming an interior initial point and for updating
% the output data in the corresponding data model element for in-service
% buses from the math model solution.

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
        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %
            vv = mm.get_idx();
            varef1 = mm.interior_va(nm, dm);
            vm = obj.interior_vm(mm, nm, dm);
            x0(vv.i1.Va:vv.iN.Va) = varef1; %% angles set to 1st ref angle
            x0(vv.i1.Vm:vv.iN.Vm) = vm;     %% voltage magnitudes
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on voltage magnitudes
            vv = mm.get_idx('var');
            lambda = mm.soln.lambda;
            mu_vm_lb = lambda.lower(vv.i1.Vm:vv.iN.Vm);
            mu_vm_ub = lambda.upper(vv.i1.Vm:vv.iN.Vm);

            %% shadow prices on node power balance
            [lam_p, lam_q] = mm.node_power_balance_prices(nm);
            lam_p = lam_p(nn.i1.bus:nn.iN.bus);     %% for bus nodes only
            lam_q = lam_q(nn.i1.bus:nn.iN.bus);     %% for bus nodes only

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
