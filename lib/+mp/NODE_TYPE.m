classdef (Sealed) NODE_TYPE
%mp.NODE_TYPE  Defines enumerated type for node types.
%
%   mp.NODE_TYPE.PQ = 1
%   mp.NODE_TYPE.PV = 2
%   mp.NODE_TYPE.REF = 3
%   mp.NODE_TYPE.NONE = 4

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties (Constant)
        PQ = 1;
        PV = 2;
        REF = 3;
        NONE = 4;
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
            TorF = val == mp.NODE_TYPE.PQ  | val == mp.NODE_TYPE.PV | ...
                   val == mp.NODE_TYPE.REF | val == mp.NODE_TYPE.NONE;
        end
    end
end
