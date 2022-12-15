classdef (Abstract) nme_gen3p < mp.nm_element % & mp.form_ac

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'gen3p';
        end

        function np = np(obj)
            np = 3;     %% this is a 3 port element
        end

        function nz = nz(obj)
            nz = 3;     %% 3 (possibly complex) non-voltage states per element
        end

        function obj = add_zvars(obj, nm, dm, idx)
            p = idx{1};
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            pg_start = dme.(sprintf('pg%d_start', p));
            qg_start = dme.(sprintf('qg%d_start', p));

            if p == 1
                nm.init_indexed_name('zr', 'Pg3', {obj.nz});
                nm.init_indexed_name('zi', 'Qg3', {obj.nz});
            end
            nm.add_var('zr', 'Pg3', {p}, ng, pg_start, 0, Inf);
            nm.add_var('zi', 'Qg3', {p}, ng, qg_start, -Inf, Inf);
        end

        function obj = build_params(obj, nm, dm)
            build_params@mp.nm_element(obj, nm, dm);    %% call parent
            obj.N = -speye(obj.nk * obj.nz);
        end
    end     %% methods
end         %% classdef
