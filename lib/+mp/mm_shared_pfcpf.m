classdef (Abstract) mm_shared_pfcpf < handle
% mp.mm_shared_pfcpf - Mixin class for PF/CPF **math model** objects.
%
% An abstract mixin class inherited by all power flow (PF) and continuation
% power flow (CPF) **math model** objects.

%   MATPOWER
%   Copyright (c) 2022-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function ad = build_aux_data(obj, nm, dm, mpopt)
            %

            %% create aux_data struct
            ad = obj.build_base_aux_data(nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
