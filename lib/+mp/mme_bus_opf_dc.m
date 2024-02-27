classdef mme_bus_opf_dc < mp.mme_bus
% mp.mme_bus_opf_dc - Math model element for bus for DC OPF.
%
% Math model element class for bus elements for DC OPF problems.
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
            x0(vv.i1.Va:vv.iN.Va) = varef1; %% angles set to 1st ref angle
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            %% bus voltage angles
            nn = nm.get_idx('node');
            va = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on node power balance
            ll = mm.get_idx('lin');
            lambda = mm.soln.lambda;
            lam_p = lambda.mu_u(ll.i1.Pmis:ll.iN.Pmis) - ...
                    lambda.mu_l(ll.i1.Pmis:ll.iN.Pmis);
            lam_p = lam_p(nn.i1.bus:nn.iN.bus);     %% for bus nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = va * 180/pi;
            dme.tab.vm(dme.on) = 1;
            dme.tab.lam_p(dme.on) = lam_p / dm.base_mva;
        end
    end     %% methods
end         %% classdef
