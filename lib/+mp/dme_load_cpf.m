classdef dme_load_cpf < mp.dme_load
%MP.DME_LOAD_CPF  MATPOWER data model class for load data

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
            vars = horzcat( export_vars@mp.dme_load(obj), ...
                {'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z'} );
        end

        function dm = parameterized(obj, dm, dmb, dmt, lam)
            load = dm.elements.load;
            b = dmb.elements.load.tab;      %% base load table
            t = dmt.elements.load.tab;      %% target load table

            load.tab.pd = b.pd + lam * (t.pd - b.pd);
            load.tab.qd = b.qd + lam * (t.qd - b.qd);
            load.tab.pd_i = b.pd_i + lam * (t.pd_i - b.pd_i);
            load.tab.qd_i = b.qd_i + lam * (t.qd_i - b.qd_i);
            load.tab.pd_z = b.pd_z + lam * (t.pd_z - b.pd_z);
            load.tab.qd_z = b.qd_z + lam * (t.qd_z - b.qd_z);
        end
    end     %% methods
end         %% classdef
