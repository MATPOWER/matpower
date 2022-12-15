classdef (Abstract) net_model_ac < mp.net_model% & mp.form_ac
%MP.NET_MODEL_AC Abstract class, explicitly a subclass of MP.NET_MODEL and
%              implicitly assumed to be a subclass of MP.FORM_AC as well

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        zr = [];
        zi = [];
    end

    properties (Access=protected)
        inln_list = {};         %% list of indexes of nme's w/inln
        snln_list = {};         %% list of indexes of nme's w/snln
        inln_hess_list = {};    %% list of indexes of nme's w/inln_hess
        snln_hess_list = {};    %% list of indexes of nme's w/snln_hess
    end

    methods
        function obj = def_set_types(obj)
            def_set_types@mp.net_model(obj);    %% call parent first
            obj.set_types.zr = 'NON-VOLTAGE VARS REAL (zr)';
            obj.set_types.zi = 'NON-VOLTAGE VARS IMAG (zi)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mp.net_model(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.Y = obj.stack_matrix_params('Y', 1);
            obj.L = obj.stack_matrix_params('L', 0);
            obj.M = obj.stack_matrix_params('M', 1);
            obj.N = obj.stack_matrix_params('N', 0);
            obj.i = obj.stack_vector_params('i');
            obj.s = obj.stack_vector_params('s');

            %% add general nonlinear function if any element has one defined
            for k = 1:length(obj.elements)
                if ~isempty(obj.elements{k}.inln)
                    obj.inln_list{end+1} = k;
                    if ~isempty(obj.elements{k}.inln_hess)
                        obj.inln_hess_list{end+1} = k;
                    end
                end
                if ~isempty(obj.elements{k}.snln)
                    obj.snln_list{end+1} = k;
                    if ~isempty(obj.elements{k}.snln_hess)
                        obj.snln_hess_list{end+1} = k;
                    end
                end
            end
            if ~isempty(obj.inln_list)
                obj.inln = @(x_, sysx, idx)port_inj_nln(obj, 'i', x_, sysx, idx);
                if ~isempty(obj.inln_hess_list)
                    obj.inln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 'i', x_, lam, sysx, idx);
                end
            end
            if ~isempty(obj.snln_list)
                obj.snln = @(x_, sysx, idx)port_inj_nln(obj, 's', x_, sysx, idx);
                if ~isempty(obj.snln_hess_list)
                    obj.snln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 's', x_, lam, sysx, idx);
                end
            end
        end

        function [g, gv1, gv2, gzr, gzi] = port_inj_nln(obj, si, x_, sysx, idx)
            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln'];
            fcn_list = [fcn '_list'];

            %% initialize
            if isempty(idx)
                sel = 0;        %% all ports
                np = obj.np;
            else
                sel = 1;        %% selected ports only
                np = length(idx);
            end
            nv = obj.get_nv_(sysx);
            nz = obj.nz;
            nc = size(x_, 2);   %% num of cols in x_, for evaluating multiple x_
            g = zeros(np, nc);
            gv1 = sparse(np, nv);
            gv2 = sparse(np, nv);
            gzr = sparse(np, nz);
            gzi = sparse(np, nz);

            %% loop through elements w/gen nonlin fcns, evaluate them
            pp = obj.get_idx('port');
            if ~sysx
                ss = obj.get_idx('state');
            end
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.elements
                nme = obj.elements{k};
                i1 = pp.i1.(nme.name)(1);
                iN = pp.iN.(nme.name)(end);

                %% set up port index vector for nme
                if sel
                    apidx = find(idx >= i1 & idx <= iN);   %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    nme_idx = idx(apidx) - i1 + 1;  %% port index vector for nme
                else
                    nme_idx = [];                   %% all ports for nme
                end

                %% set up proper x_ for nme
                if sysx
                    nme_x_ = x_;
                else
                    if isfield(ss.i1, nme.name)
                        j1 = ss.i1.(nme.name)(1);
                        jN = ss.iN.(nme.name)(end);
                    else
                        j1 = 1;
                        jN = 0;
                    end
                    nme_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];
                end

                %% call nonlinear function
                gg = cell(1, nargout);
                [gg{:}] = nme.(fcn)(nme_x_, sysx, nme_idx);

                %% insert the results in aggregate output args
                if sel
                    g(apidx, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(apidx, :) = gg{2};
                            gv2(apidx, :) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(apidx, :) = gg{4};
                                gzi(apidx, :) = gg{5};
                            end
                        else
                            gv1(apidx, i1:iN) = gg{2};
                            gv2(apidx, i1:iN) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(apidx, j1:jN) = gg{4};
                                gzi(apidx, j1:jN) = gg{5};
                            end
                        end
                    end
                else
                    g(i1:iN, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(i1:iN, :) = gg{2};
                            gv2(i1:iN, :) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(i1:iN, :) = gg{4};
                                gzi(i1:iN, :) = gg{5};
                            end
                        else
                            gv1(i1:iN, i1:iN) = gg{2};
                            gv2(i1:iN, i1:iN) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(i1:iN, j1:jN) = gg{4};
                                gzi(i1:iN, j1:jN) = gg{5};
                            end
                        end
                    end
                end
            end     %% for loop
        end

        function H = port_inj_nln_hess(obj, si, x_, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x_, lam)
            % H = obj.port_inj_power_hess(x_, lam, sysx)
            % H = obj.port_inj_power_hess(x_, lam, sysx, idx)
            if nargin < 6
                idx = [];
                if nargin < 5
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln_hess'];
            fcn_list = [fcn '_list'];

            %% initialize
            n = 2 * length(x_);
            H = sparse(n, n);

            %% loop through elements w/gen nonlin Hessians, evaluate them
            pp = obj.get_idx('port');
            if ~sysx
                ss = obj.get_idx('state');
            end
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.elements
                nme = obj.elements{k};
                i1 = pp.i1.(nme.name)(1);
                iN = pp.iN.(nme.name)(end);

                %% set up x_ for nme & corresp row/col indices for nme
                if sysx
                    nme_x_ = x_;
                else
                    nv = obj.get_nv_(sysx);
                    nz = obj.nz;
                    if isfield(ss.i1, nme.name)
                        j1 = ss.i1.(nme.name)(1);
                        jN = ss.iN.(nme.name)(end);
                    else
                        j1 = 1;
                        jN = 0;
                    end
                    nme_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];

                    %% indices of rows/cols of H corresponding to nme x_
                    h = [(i1:iN) nv+(i1:iN) 2*nv+(j1:jN) 2*nv+nz+(j1:jN)].';
                end

                %% set up port index and lambda vectors for nme
                if ~isempty(idx)    %% selected ports only
                    apidx = find(idx >= i1 & idx <= iN);    %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    nme_idx = idx(apidx) - i1 + 1;  %% port index vector for nme
                    nme_lam = lam(apidx);           %% corresponding lam
                else                %% all ports
                    nme_idx = [];
                    nme_lam = lam(i1:iN);
                end

                %% call nonlinear function
                nme_H = nme.(fcn)(nme_x_, nme_lam, sysx, nme_idx);

                %% accumulate output
                if sysx
                    H = H + nme_H;
                else
                    H(h,h) = H(h,h) + nme_H;
                end
            end
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_current_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1);
                Gv1 = C * Iv1;      %% Gva or Gvr
                Gv2 = C * Iv2;      %% Gvm or Gvi
                Gzr = C * Izr;
                Gzi = C * Izi;
            else
                I = obj.port_inj_current(x_, 1);
            end

            %% nodal current balance
            G = C * I;
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_power_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1);
                Gv1 = C * Sv1;      %% Gva or Gvr
                Gv2 = C * Sv2;      %% Gvm or Gvi
                Gzr = C * Szr;
                Gzi = C * Szi;
            else
                S = obj.port_inj_power(x_, 1);
            end

            %% nodal power balance
            G = C * S;
        end

        function d2G = nodal_complex_current_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_current_hess(x_, obj.C' * lam);
        end

        function d2G = nodal_complex_power_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x_, obj.C' * lam);
        end

        function obj = port_inj_soln(obj)
            %% compute port injections
%             obj.soln.gi_ = obj.port_inj_current(obj.soln.x);
            obj.soln.gs_ = obj.port_inj_power(obj.soln.x);
        end

        function va = get_va(obj, idx)
            if isfield(obj.soln, 'v')           %% solved value
                if nargin < 2 || isempty(idx)
                    va = angle(obj.soln.v);
                else
                    va = angle(obj.soln.v(idx));
                end
            else                                %% initial value
                if nargin < 2 || isempty(idx)
                    va = obj.initial_voltage_angle();
                else
                    va = obj.initial_voltage_angle(idx);
                end
            end
        end
    end     %% methods
end         %% classdef
