classdef nme_bus_acc < mp.nme_bus & mp.form_acc
% mp.nme_bus_acc - Network model element for bus for AC cartesian voltage formulations.
%
% Adds voltage variables ``Vr`` and ``Vi`` to the network model and inherits
% from mp.form_acc.

%   MATPOWER
%   Copyright (c) 2018-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
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

            %% prepare box bounds for voltage coordinates
            V0 = dme.vm_start .* exp(1j * dme.va_start);
            vclim = 1.1 * dme.vm_ub;

            nm.add_var('vr', 'Vr', nb, real(V0), -vclim, vclim);
            nm.add_var('vi', 'Vi', nb, imag(V0), -vclim, vclim);
        end
    end     %% methods
end         %% classdef
