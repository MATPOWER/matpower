classdef mme_reserve_zone < mp.mm_element
% mp.mme_reserve_zone - Mathematical model element for reserve zone.
%
% Math model element class for reserve zone elements.
%
% Implements methods for adding reserve zone constraints, and for updating
% the output data in the corresponding data model element for in-service
% reserve zones from the math model solution.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            %
            name = 'reserve_zone';
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            dme = obj.data_model_element(dm);
            rgen_dme = dm.elements.reserve_gen;
            igr = rgen_dme.gen; %% index of online gen for each online reserve gen

            Areq = dme.zones(:, igr);
            lreq = dme.req;

            mm.add_lin_constraint('Rreq', Areq, lreq, [], {'R'});
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            %

            dme = obj.data_model_element(dm);

            %% get reserve zone prices
            ll = mm.get_idx('lin');
            prc = mm.soln.lambda.mu_l(ll.i1.Rreq:ll.iN.Rreq);

            %% update in the data model
            dme.tab.prc(dme.on) = prc / dm.base_mva;
        end
    end     %% methods
end         %% classdef
