classdef (Abstract) mm_shared_pfcpf_ac_i < handle
% mp.mm_shared_pfcpf_ac_i - Mixin class for AC-current PF/CPF **math model** objects.
%
% An abstract mixin class inherited by all AC power flow (PF) and continuation
% power flow (CPF) **math model** objects that use a current balance
% formuation.
%
% Code shared between AC cartesian and polar formulations with current balance
% belongs in this class.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function ad = build_aux_data_i(obj, nm, ad)
            %

            %% build additional aux data
            N = nm.C(ad.pv, :) * nm.N;      %% z coefficients for z @ PV nodes
            [ii, jj, ~] = find(N == -1);    %% find element representing 1st
            [~, ia] = unique(ii, 'first');  %% direct injection at each PV node
            [~, ib] = sort(ii(ia));         %% sort by PV node

            %% save additional aux data
            ad.zi_idx = jj(ia(ib));         %% z-var index corresp to PV nodes
        end
    end     %% methods
end         %% classdef
