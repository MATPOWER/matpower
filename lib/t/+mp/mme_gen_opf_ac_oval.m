classdef mme_gen_opf_ac_oval < mp.mme_gen_opf_ac
% mp.mme_gen_opf_ac_oval - Math model element for generator for AC OPF w/oval cap curve.
%
% Math model element class for generator elements for AC OPF problems,
% implementing an oval, as opposed to rectangular, PQ capability curve.

%   MATPOWER
%   Copyright (c) 2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            % Set up the nonlinear constraint for gen oval PQ capability curves.
            % ::
            %
            %   mme.add_constraints(mm, nm, dm, mpopt)

            dme = obj.data_model_element(dm);

            %% generator PQ capability curve constraints
            idx = [];       %% which generators get this constraint
                            %% empty ==> all
            if isempty(idx)
                idx = (1:dme.n)';
            end

            %% get generator limit data
            p_lb = dme.pg_lb(idx);
            p_ub = dme.pg_ub(idx);
            q_lb = dme.qg_lb(idx);
            q_ub = dme.qg_ub(idx);

            %% compute oval specs, all vectors, 4 params per gen
            a2 = (p_ub - p_lb) .^ 2;        % square of horizontal (p) radius
            b2 = ((q_ub - q_lb) / 2) .^ 2;  % square of vertical (q) radius
            p0 = p_lb;                      % horizontal (p) center
            q0 = (q_ub + q_lb) / 2;         % vertical (q) center

            %% add constraint
            fcn = @(xx)oval_pq_capability_fcn(obj, xx, idx, p0, q0, a2, b2);
            hess = @(xx, lam)oval_pq_capability_hess(obj, xx, lam, idx, p0, q0, a2, b2);
            mm.add_nln_constraint('PQoval', dme.n, 0, fcn, hess, {'Pg', 'Qg'});

            %% call parent
            add_constraints@mp.mme_gen_opf_ac(obj, mm, nm, dm, mpopt);
        end

        function [h, dh] = oval_pq_capability_fcn(obj, xx, idx, p0, q0, a2, b2)
            % Compute oval PQ capability constraints and Jacobian.
            % ::
            %
            %   h = mme.oval_pq_capability_fcn(xx, idx, p0, q0, a2, b2)
            %   [h, dh] = mme.oval_pq_capability_fcn(xx, idx, p0, q0, a2, b2)
            %
            % Compute constraint function and optionally the Jacobian for
            % oval PQ capability limits.
            %
            % Inputs:
            %   xx (1 x 2 cell array) : active power injection in
            %       ``xx{1}``, reactive injection in ``xx{2}``
            %   idx (integer) : index of subset of generators of interest to
            %       include in constraint; if empty, include all
            %   p0 (double) : vector of horizontal (p) centers
            %   q0 (double) : vector of vertical (q) centers
            %   a2 (double) : vector of squares of horizontal (p) radii
            %   b2 (double) : vector of squares of vertical (q) radii
            %
            % Outputs:
            %   h (double) : constraint function, :math:`\h(\x)`
            %   dh (double) : constraint Jacobian, :math:`\h_\x`
            %
            % Note that the oval specs ``p0``, ``q0``, ``a2``, ``b2``
            % are assumed to have dimension corresponding to ``idx``.

            [p, q] = deal(xx{:});
            ng = length(p);
            if ~isempty(idx)
                p = p(idx);
                q = q(idx);
            end

            %% evaluate constraint function
            h = (p - p0).^2 ./ a2 + (q - q0).^2 ./ b2 - 1;

            %% evaluate constraint Jacobian
            if nargout > 1
                dhdp = spdiags(2*(p - p0) ./ a2, 0, ng, ng);
                dhdq = spdiags(2*(q - q0) ./ b2, 0, ng, ng);
                dh = [dhdp dhdq];
            end
        end

        function d2H = oval_pq_capability_hess(obj, xx, lam, idx, p0, q0, a2, b2)
            % Compute oval PQ capability constraint Hessian.
            % ::
            %
            %   d2H = mme.oval_pq_capability_hess(xx, lam, idx, p0, q0, a2, b2)
            %
            % Compute a sparse Hessian matrix for oval PQ capability limits.
            % Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the constraint Jacobian by a vector :math:`\muv`.
            %
            % Inputs:
            %   xx (1 x 2 cell array) : active power injection in
            %       ``xx{1}``, reactive injection in ``xx{2}``
            %   lam (double) : vector :math:`\muv` of multipliers
            %   idx (integer) : index of subset of generators of interest to
            %       include in constraint; if empty, include all
            %   p0 (double) : vector of horizontal (p) centers
            %   q0 (double) : vector of vertical (q) centers
            %   a2 (double) : vector of squares of horizontal (p) radii
            %   b2 (double) : vector of squares of vertical (q) radii
            %
            % Output:
            %   d2H (double) : sparse constraint Hessian matrix
            %
            % Note that the oval specs ``p0``, ``q0``, ``a2``, ``b2``
            % are assumed to have dimension corresponding to ``idx``.

            [p, q] = deal(xx{:});
            if ~isempty(idx)
                p = p(idx);
                q = q(idx);
            end
            ng = length(p);
            zz = sparse(ng, ng);

            %% evaluate constraint Hessian
            d2H_pp = sparse(1:ng, 1:ng, 2 * lam ./ a2, ng, ng);
            d2H_qq = sparse(1:ng, 1:ng, 2 * lam ./ b2, ng, ng);
            d2H = [ d2H_pp  zz;
                    zz      d2H_qq ];
        end
    end     %% methods
end         %% classdef
