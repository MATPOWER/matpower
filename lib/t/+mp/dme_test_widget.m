classdef dme_test_widget < mp.dm_element
%MP.DME_TEST_WIDGET  MATPOWER data model class for test widget data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        xtab
    end     %% properties

    methods
        function name = name(obj)
            name = 'test_widget';
        end

        function label = label(obj)
            label = 'Test Widget';
        end

        function label = labels(obj)
            label = 'Test Widgets';
        end

        function names = main_table_var_names(obj)
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'ids', 'twos', 'alpha', 'beta', 'gamma', 'delta', 'epsilon'});
        end

        function vars = export_vars(obj)
            vars = {'alpha', 'gamma', 'epsilon'};
        end

        function s = export_vars_offline_val(obj)
            %

            s = export_vars_offline_val@mp.dm_element(obj);     %% call parent
            s.alpha = 2;
            s.epsilon = 0;
        end
    end     %% methods
end         %% classdef
