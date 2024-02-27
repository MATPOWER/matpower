classdef nme_bus3p_acc < mp.nme_bus3p & mp.form_acc
% mp.nme_bus3p_acc - Network model element for 3-phase bus, AC cartesian voltage formulation.
%
% Adds voltage variables ``Vr3`` and ``Vi3`` to the network model and inherits
% from mp.form_acc.

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

            %% prepare box bounds for voltage coordinates
            vm_start = dme.(sprintf('vm%d_start', p));
            va_start = dme.(sprintf('va%d_start', p));
            v_start = vm_start .* exp(1j * va_start);
            vclim = 1.5;

            if p == 1
                nm.init_indexed_name('vr', 'Vr3', {obj.nn});
                nm.init_indexed_name('vi', 'Vi3', {obj.nn});
            end
            nm.add_var('vr', 'Vr3', {p}, nb, real(v_start), -vclim, vclim);
            nm.add_var('vi', 'Vi3', {p}, nb, imag(v_start), -vclim, vclim);
        end
    end     %% methods
end         %% classdef
