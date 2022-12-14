classdef nme_gizmo_acc_nln < mp.nme_gizmo_acc & mp.nme_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nme_class = @mp.nme_gizmo_acc;
    end

    methods
        function obj = nme_gizmo_acc_nln()
            obj@mp.nme_gizmo_acc();
            obj.nme_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@mp.nme_gizmo_acc(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.nme.count(dm);
            nk = count@mp.nme_gizmo_acc(obj, dm);
        end
    end     %% methods
end         %% classdef
