classdef (Abstract) nme_gen_ac < mp.nme_gen% & mp.form_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('zr', 'Pg', ng, dme.pg_start, dme.pg_lb, dme.pg_ub);
            nm.add_var('zi', 'Qg', ng, dme.qg_start, dme.qg_lb, dme.qg_ub);
        end

        function obj = build_params(obj, nm, dm)
            build_params@mp.nme_gen(obj, nm, dm);   %% call parent
            obj.N = -speye(obj.nk * obj.nz);
        end
    end     %% methods
end         %% classdef
