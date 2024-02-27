classdef math_model_cpf_acps < mp.math_model_cpf_acp & mp.mm_shared_pfcpf_acps
% mp.math_model_cpf_acps - CPF **math model** for AC-polar-power formulation.
%
% Implements formulation-specific and CPF-specific node balance constraint.
%
% Provides methods for warm-starting solver with updated data.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function tag = form_tag(obj)
            %

            tag = 'acps';
        end

        function name = form_name(obj)
            %

            name = 'AC-polar-power';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %

            %% power balance constraints
            ad = obj.aux_data;
            fcn = @(x)node_balance_equations_cpf(obj, x, nm);
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function varargout = expand_z_warmstart(obj, nm, ad, varargin)
            %

            %% expand input tangent z vectors to all nodes + lambda
            varargout = cell(size(varargin));
            i = [ad.pv; ad.pq; nm.nv/2 + ad.pq; nm.nv+1];
            for k = 1:length(varargin)
                z = zeros(nm.nv, 1);
                z(i) = varargin{k};
                varargout{k} = z;
            end
        end

        function opt = solve_opts_warmstart(obj, opt, ws, nm)
            %

            ad = obj.aux_data;

            %% update warm start states and tangent vectors
            ws.x  = [angle(ws.cV([ad.pv; ad.pq])); abs(ws.cV(ad.pq)); ws.clam];
            ws.xp = [angle(ws.pV([ad.pv; ad.pq])); abs(ws.pV(ad.pq)); ws.plam];
            opt.x0 = ws.x;   %% ignored, overridden by ws.x

            %% reduce tangent vectors for this mm
            i = [ad.pv; ad.pq; nm.nv/2 + ad.pq; nm.nv+1];
            ws.z  = ws.z(i);
            ws.zp = ws.zp(i);
            opt.warmstart = ws;
        end
    end     %% methods
end         %% classdef
