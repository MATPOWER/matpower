classdef (Sealed) NODE_TYPE
%mp.NODE_TYPE - Defines enumerated type for node types.
%
% mp.NODE_TYPE Properties:
%   * PQ - PQ node (= 1)
%   * PV - PV node (= 2)
%   * REF - reference node (= 3)
%   * NONE - isolated node (= 4)
%
% mp.NODE_TYPE Methods:
%   * is_valid - returns true if the value is a valid node type
%
% All properties are ``Constant`` properties and the class is a ``Sealed``
% class. So the properties function as global constants which do not create
% an instance of the class, e.g. mp.NODE_TYPE.REF.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties (Constant)
        PQ = 1;     % PQ node
        PV = 2;     % PV node
        REF = 3;    % reference node
        NONE = 4;   % isolated node
    end

%     methods (Access = private)      %% to prevent instantiation
%         function obj = NODE_TYPE
%         end
%     end

    methods (Static)
%         function TorF = is_pq(val)
%             TorF = val == mp.NODE_TYPE.PQ;
%         end
% 
%         function TorF = is_pv(val)
%             TorF = val == mp.NODE_TYPE.PV;
%         end
% 
%         function TorF = is_ref(val)
%             TorF = val == mp.NODE_TYPE.REF;
%         end

        function TorF = is_valid(val)
            % Returns true if the value is a valid node type.
            % ::
            %
            %   TorF = mp.NODE_TYPE.is_valid(val)
            %
            % Input:
            %   val (integer) : node type value to check for validity
            %
            % Output:
            %   TorF (boolean) : true if ``val`` is a valid node type

            TorF = val == mp.NODE_TYPE.PQ  | val == mp.NODE_TYPE.PV | ...
                   val == mp.NODE_TYPE.REF | val == mp.NODE_TYPE.NONE;
        end
    end
end
