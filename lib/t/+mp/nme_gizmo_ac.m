classdef (Abstract) nme_gizmo_ac < mp.nme_gizmo% & mp.form_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = add_zvars(obj, nm, dm, idx)
            tab = obj.data_model_element(dm).tab;
            nk = obj.nk;
            switch idx{:}
                case 1
                    Zmax = ones(nk, 1);
                    Zr   = tab.Zr1;
                    Zi   = tab.Zi1;
                case 2
                    Zmax = 2 * ones(nk, 1);
                    Zr   = tab.Zr2;
                    Zi   = tab.Zi2;
            end
            vname_r = sprintf('Zr%d_gizmo', idx{:});
            vname_i = sprintf('Zi%d_gizmo', idx{:});
            nm.add_var('zr', vname_r, nk, Zr, -Zmax, Zmax);
            nm.add_var('zi', vname_i, nk, Zi, -Zmax, Zmax);
        end

        function obj = build_params(obj, nm, dm)
            build_params@mp.nme_gizmo(obj, nm, dm);    %% call parent
            tab = obj.data_model_element(dm).tab;
            nk = obj.nk;

            %% collect parameters from data table
            y1 = tab.Y1r + 1j * tab.Y1i;
            y2 = tab.Y2r + 1j * tab.Y2i;
            ll = tab.Lr + 1j * tab.Li;
            ii = tab.Ir + 1j * tab.Ii;
            m1 = tab.M1r + 1j * tab.M1i;
            m2 = tab.M2r + 1j * tab.M2i;
            nn = tab.Nr + 1j * tab.Ni;
            ss = tab.Sr + 1j * tab.Si;
            zz = zeros(nk, 1);

            %% construct model parameters
            j1 = (1:nk);
            j2 = nk+j1;
            j3 = nk+j2;
            obj.Y = sparse( ...
                [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                [y1; zz; -y1; zz; y2; zz; -y1; zz; y1], 3*nk, 3*nk );
            obj.L = sparse( ...
                [j1 j1 j2 j2 j3 j3 ]', ...
                [j1 j2 j1 j2 j1 j2 ]', ...
                [zz; ll; zz; -ll; zz; zz], 3*nk, 2*nk );
            obj.i = [-ii; ii; zz];
            obj.M = sparse( ...
                [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                [m1; -m1; zz; -m1; m1; zz; zz; zz; m2], 3*nk, 3*nk );
            obj.N = sparse( ...
                [j1 j1 j2 j2 j3 j3 ]', ...
                [j1 j2 j1 j2 j1 j2 ]', ...
                [zz; zz; nn; zz; -nn; zz], 3*nk, 2*nk );
            obj.s = [zz; -ss; ss];
        end
    end     %% methods
end         %% classdef
