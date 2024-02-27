classdef mme_bus3p_opf_acc < mp.mme_bus3p
% mp.mme_bus3p_opf_acc - Math model element for 3-phase bus for AC cartesian voltage OPF.
%
% Math model element class for 3-phase bus elements for AC cartesian voltage
% OPF problems.
%
% Implements method for forming an interior initial point.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
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
            vm = 1;     %% voltage magnitude set to 1 p.u.
            %% voltage angles set to angle of 1st ref node
            %% + phase offset, i.e. +/-120 deg
            varef1 = mm.interior_va(nm, dm);
            v1 = vm * exp(1j*varef1);
            v2 = vm * exp(1j*(varef1-2*pi/3));
            v3 = vm * exp(1j*(varef1+2*pi/3));

            x0(vv.i1.Vr3(1):vv.iN.Vr3(1)) = real(v1);
            x0(vv.i1.Vr3(2):vv.iN.Vr3(2)) = real(v2);
            x0(vv.i1.Vr3(3):vv.iN.Vr3(3)) = real(v3);
            x0(vv.i1.Vi3(1):vv.iN.Vi3(1)) = imag(v1);
            x0(vv.i1.Vi3(2):vv.iN.Vi3(2)) = imag(v2);
            x0(vv.i1.Vi3(3):vv.iN.Vi3(3)) = imag(v3);
        end
    end     %% methods
end         %% classdef
