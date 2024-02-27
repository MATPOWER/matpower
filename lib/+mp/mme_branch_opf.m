classdef (Abstract) mme_branch_opf < mp.mme_branch
% mp.mme_branch_opf - Math model element abstract base class for branch for OPF.
%
% Math model element abstract base class for branch elements, including
% transmission lines and transformers, for OPF problems.
%
% Implements methods to prepare data required for angle difference limit
% constraints and to extract shadow prices for these constraints from
% the math model solution.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function [A, l, u, i] = ang_diff_params(obj, dm, ignore)
            %
            dme = obj.data_model_element(dm);

            %% from makeAang()
            nb = dm.elements.bus.n;
            branch = dme.tab;

            if ignore
                A  = sparse(0, nb);
                l  = [];
                u  = [];
                i  = [];
            else
                i = find( ...
                    branch.vad_lb ~= 0 & ...
                        (branch.vad_lb > -360 | branch.vad_ub == 0) | ...
                    branch.vad_ub ~= 0 & ...
                        (branch.vad_ub <  360 | branch.vad_lb == 0) );
                n = length(i);

                if n > 0
                    ii = [(1:n)'; (1:n)'];
                    jj = [dme.fbus(i); dme.tbus(i)];
                    A = sparse(ii, jj, [ones(n, 1); -ones(n, 1)], n, nb);
                    l = branch.vad_lb(i);
                    u = branch.vad_ub(i);
                    l(l < -360) = -Inf;
                    u(u >  360) =  Inf;
                    l = l * pi/180;
                    u = u * pi/180;
                else
                    A = sparse(0, nb);
                    l =[];
                    u =[];
                end
            end
        end

        function [mu_vad_lb, mu_vad_ub] = ang_diff_prices(obj, mm, nme)
            %

            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(nme.nk, 1);
            mu_vad_ub = mu_vad_lb;
            if length(iang)
                ll = mm.get_idx('lin');
                lambda = mm.soln.lambda;
                mu_vad_lb(iang) = lambda.mu_l(ll.i1.ang:ll.iN.ang);
                mu_vad_ub(iang) = lambda.mu_u(ll.i1.ang:ll.iN.ang);
            end
        end
    end     %% methods
end         %% classdef
