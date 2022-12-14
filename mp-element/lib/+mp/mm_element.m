classdef (Abstract) mm_element < handle
%MP.MM_ELEMENT  Abstract base class for mathematical model element

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function nme = network_model_element(obj, nm, name)
            if nargin < 3
                name = obj.name;
            end
            nme = nm.elements.(name);
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
        end
    end     %% methods
end         %% classdef
