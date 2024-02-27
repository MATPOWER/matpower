classdef (Abstract) nme_branch_ac < mp.nme_branch% & mp.form_ac
% mp.nme_branch_ac - Network model element abstract base class for branch for AC formulations.
%
% Implements building of the admittance parameter :math:`\YY` for branches.

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
            % Builds the admittance parameter :math:`\YY` for branches.

            build_params@mp.nme_branch(obj, nm, dm);    %% call parent

            dme = obj.data_model_element(dm);
            nl = obj.nk;

            tm = ones(nl, 1);           %% default tap ratio = 1
            i = find(dme.tm);           %% indices of non-zero tap ratios
            tm(i) = dme.tm(i);              %% assign non-zero tap ratios
            T = tm .* exp(1j * dme.ta);     %% add phase shifters

            ys = 1 ./ (dme.r + 1j * dme.x);     %% series admittance
            yff = (ys + dme.g_fr + 1j * dme.b_fr) ./ (T .* conj(T));
            ytt = ys + dme.g_to + 1j * dme.b_to;
            yft = - ys ./ conj(T);
            ytf = - ys ./ T;

            obj.Y = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [yff; yft; ytf; ytt], 2*nl, 2*nl );

% same as:
%             Yff = spdiags(yff, 0, nl, nl);
%             Yft = spdiags(yft, 0, nl, nl);
%             Ytf = spdiags(ytf, 0, nl, nl);
%             Ytt = spdiags(ytt, 0, nl, nl);
%             obj.Y = [ Yff Yft;
%                       Ytf Ytt  ];
        end
    end     %% methods
end         %% classdef
