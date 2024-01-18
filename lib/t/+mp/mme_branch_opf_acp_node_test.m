classdef mme_branch_opf_acp_node_test < mp.mme_branch_opf_acp

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function [A, l, u, i] = ang_diff_params(obj, dm, ignore)
            dme = obj.data_model_element(dm);

            %% from makeAang()
            nb = 0;
            for k = 1:dme.nbet
                if dm.elements.has_name(dme.cxn_type{k})
                    nb = nb + dm.elements.(dme.cxn_type{k}).n;
                end
            end
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
            if length(i)
                warning('OPF branch angle difference limits not implemented for this case.');
            end
        end
    end     %% methods
end         %% classdef
