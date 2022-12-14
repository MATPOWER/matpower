classdef dme_shunt_cpf < mp.dme_shunt
%MP.DME_SHUNT_CPF  MATPOWER data model class for shunt data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function vars = export_vars(obj)
            vars = horzcat( export_vars@mp.dme_shunt(obj), {'gs', 'bs'} );
        end

        function dm = parameterized(obj, dm, dmb, dmt, lam)
            shunt = dm.elements.shunt;
            b = dmb.elements.shunt.tab;     %% base shunt table
            t = dmt.elements.shunt.tab;     %% target shunt table

            shunt.tab.gs = b.gs + lam * (t.gs - b.gs);
            shunt.tab.bs = b.bs + lam * (t.bs - b.bs);
        end
    end     %% methods
end         %% classdef
