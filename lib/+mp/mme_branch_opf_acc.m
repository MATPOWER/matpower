classdef mme_branch_opf_acc < mp.mme_branch_opf_ac
% mp.mme_branch_opf_acc - Math model element for branch for AC cartesian voltage OPF.
%
% Math model element class for branch elements, including transmission lines
% and transformers, for AC cartesian voltage OPF problems.
%
% Implements method for adding branch angle difference constraints and
% overrides method to extract shadow prices for these constraints from
% the math model solution.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            nme = obj.network_model_element(nm);

            %% call parent
            add_constraints@mp.mme_branch_opf_ac(obj, mm, nm, dm, mpopt);

            %% branch angle difference limits
            [Aang, lang, uang, iang] = obj.ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            nang = length(iang);
            if nang
                fcn_ang = @(xx)ang_diff_fcn(nme, xx, Aang, lang, uang);
                hess_ang = @(xx, lam)ang_diff_hess(nme, xx, lam, Aang);
                mm.add_nln_constraint({'angL', 'angU'}, [nang;nang], 0, fcn_ang, hess_ang, {'Vr', 'Vi'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end

        function [mu_vad_lb, mu_vad_ub] = ang_diff_prices(obj, mm, nme)
            %

            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(nme.nk, 1);
            mu_vad_ub = mu_vad_lb;
            if length(iang)
                nni = mm.get_idx('nli');
                lambda = mm.soln.lambda;
                mu_vad_ub(iang) = lambda.ineqnonlin(nni.i1.angU:nni.iN.angU);
                mu_vad_lb(iang) = lambda.ineqnonlin(nni.i1.angL:nni.iN.angL);
            end
        end
    end     %% methods
end         %% classdef
