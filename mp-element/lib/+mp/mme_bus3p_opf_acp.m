classdef mme_bus3p_opf_acp < mp.mme_bus3p

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
        function x0 = interior_x0(obj, mm, nm, dm, x0)
            vv = mm.get_idx();
            vm = 1;     %% voltage magnitude set to 1 p.u.
            %% voltage angles set to angle of 1st ref node
            %% + phase offset, i.e. +/-120 deg
            varef1 = mm.interior_va(nm, dm);

            x0(vv.i1.Va3(1):vv.iN.Va3(1)) = varef1;
            x0(vv.i1.Va3(2):vv.iN.Va3(2)) = varef1-2*pi/3;
            x0(vv.i1.Va3(3):vv.iN.Va3(3)) = varef1+2*pi/3;
            x0(vv.i1.Vm3(1):vv.iN.Vm3(1)) = vm;
            x0(vv.i1.Vm3(2):vv.iN.Vm3(2)) = vm;
            x0(vv.i1.Vm3(3):vv.iN.Vm3(3)) = vm;
        end
    end     %% methods
end         %% classdef
