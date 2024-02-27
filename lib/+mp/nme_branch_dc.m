classdef nme_branch_dc < mp.nme_branch & mp.form_dc
% mp.nme_branch_dc - Network model element for branch for DC formulations.
%
% Implements building of the branch parameters :math:`\BB` and
% :math:`\pv`, and inherits from :class:`mp.form_dc`.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_branch(obj, nm, dm);    %% call parent

            dme = obj.data_model_element(dm);
            nl = obj.nk;

            tm = ones(nl, 1);           %% default tap ratio = 1
            i = find(dme.tm);           %% indices of non-zero tap ratios
            tm(i) = dme.tm(i);          %% assign non-zero tap ratios

            b = 1 ./ dme.x;             %% series susceptance
            b = b ./ tm;
            Pfinj = b .* (-dme.ta);
            obj.B = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [b; -b; -b; b], ...
                2*nl, 2*nl );
            obj.p = [Pfinj + dme.g_fr; -Pfinj + dme.g_to];
        end
    end     %% methods
end         %% classdef
