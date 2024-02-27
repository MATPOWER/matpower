classdef mme_bus_opf_acc < mp.mme_bus_opf_ac
% mp.mme_bus_opf_acc - Math model element for bus for AC cartesian voltage OPF.
%
% Math model element class for bus elements for AC cartesian voltage OPF
% problems.
%
% Implements methods for adding constraints for reference voltage angle,
% fixed voltage magnitudes and voltage magnitude limits, for forming an
% interior initial point and for updating the output data in the
% corresponding data model element for in-service buses from the math
% model solution.

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
        function add_constraints(obj, mm, nm, dm, mpopt)
            %
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% voltage angle reference constraint
            ref = find(nme.node_types(nm, dm) == mp.NODE_TYPE.REF);
            varef = dme.va_start(ref);
            fcn_vref = @(xx)va_fcn(nme, xx, ref, varef);
            hess_vref = @(xx, lam)va_hess(nme, xx, lam, ref);
            mm.add_nln_constraint('Vref', length(ref), 1, fcn_vref, hess_vref, {'Vr', 'Vi'});

            %% fixed voltage magnitudes
            veq = find(dme.vm_lb == dme.vm_ub);
            nveq = length(veq);
            if nveq
                fcn_vm2eq = @(xx)vm2_fcn(nme, xx, veq, dme.vm_ub(veq).^2);
                hess_vm2eq = @(xx, lam)vm2_hess(nme, xx, lam, veq);
                mm.add_nln_constraint('Veq', nveq, 1, fcn_vm2eq, hess_vm2eq, {'Vr', 'Vi'});
            end
            mm.userdata.veq = veq;

            %% voltage magnitude limits
            viq = find(dme.vm_lb ~= dme.vm_ub);
            nviq = length(viq);
            if nviq
                fcn_vlim = @(xx)vm2_fcn(nme, xx, viq, ...
                        {dme.vm_lb(viq).^2, dme.vm_ub(viq).^2} );
                hess_vlim = @(xx, lam)vm2_hess(nme, xx, lam, viq);
                mm.add_nln_constraint({'Vmin', 'Vmax'}, [nviq;nviq], 0, fcn_vlim, hess_vlim, {'Vr', 'Vi'});
            end
            mm.userdata.viq = viq;
        end

        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %
            vv = mm.get_idx();
            varef1 = mm.interior_va(nm, dm);
            vm = obj.interior_vm(mm, nm, dm);
            v_ = vm * exp(1j*varef1);
            x0(vv.i1.Vr:vv.iN.Vr) = real(v_);
            x0(vv.i1.Vi:vv.iN.Vi) = imag(v_);
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on voltage magnitudes
            [nne, nni] = mm.get_idx('nle', 'nli');
            lambda = mm.soln.lambda;
            mu_vm_lb = zeros(nn.N.bus, 1);  %% init to all 0
            mu_vm_ub = mu_vm_lb;            %% init to all 0
            if mm.userdata.veq
                lam = lambda.eqnonlin(nne.i1.Veq:nne.iN.Veq);
                lam_p = zeros(size(lam));
                lam_n = zeros(size(lam));
                lam_p(lam > 0) =  lam(lam > 0);
                lam_n(lam < 0) = -lam(lam < 0);
                mu_vm_lb(mm.userdata.veq) = lam_n;
                mu_vm_ub(mm.userdata.veq) = lam_p;
            end
            mu_vm_lb(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmin:nni.iN.Vmin);
            mu_vm_ub(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmax:nni.iN.Vmax);

            vm = abs(V);
            mu_vm_lb = mu_vm_lb .* vm * 2;
            mu_vm_ub = mu_vm_ub .* vm * 2;

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
