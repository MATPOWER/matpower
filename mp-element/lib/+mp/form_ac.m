classdef (Abstract) form_ac < mp.form
%MP.FORM_AC  MATPOWER Formulation base class for AC formulations
%   Each concrete Network Model Element class must inherit, at least
%   indirectly, from both MP.NM_ELEMENT and MP.FORM.
%
%   Subclass of MP.FORM.
%   MP.FORM provides properties and methods related to the specific
%   formulation (e.g. DC version, AC polar power version, etc.)
%
%   MP.FORM_AC defines:
%       linear current injection       = Y v_ + L z_ + i
%       linear complex power injection = M v_ + N z_ + s
%
%   Properties
%       (model parameters)
%       params - cell array of model parameter field names
%       Y - np*nk x nn matrix
%       L - np*nk x nz matrix
%       M - np*nk x nn matrix
%       N - np*nk x nz matrix
%       i - np*nk x 1 matrix
%       s - np*nk x 1 matrix
%
%   Methods
%       model_params() - cell array of names of model parameters
%                        {'Y', 'L', 'M', 'N', 'i', 's'}

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        %% model parameters
        Y = [];     %% ilin(x_, idx) = Y(idx, :) * v_ + L(idx, :) * z_ + i(idx)
        L = [];
        M = [];     %% slin(x_, idx) = M(idx, :) * v_ + N(idx, :) * z_ + s(idx)
        N = [];
        i = [];
        s = [];
        param_ncols = struct('Y', 2, 'L', 3, 'M', 2, 'N', 3, 'i', 1, 's', 1);
            %% num of columns for each parameter, where 1 = 1, 2 = np, 3 = nz
        inln = '';      %% fcn handle, i = ilin(x_, idx) + inln(x_, idx)
        snln = '';      %% fcn handle, s = slin(x_, idx) + snln(x_, idx)
        inln_hess = ''; %% fcn handle, if empty, Hessian assumed zero
        snln_hess = ''; %% fcn handle, if empty, Hessian assumed zero
    end

    methods
        function name = form_name(obj)
            name = 'AC';
        end

        function tag = form_tag(obj)
            tag = 'ac';
        end

        function params = model_params(obj)
           params = {'Y', 'L', 'M', 'N', 'i', 's'};
        end

        function vtypes = model_zvars(obj)
            vtypes = {'zr', 'zi'};
        end

        function [I, Iv1, Iv2, Izr, Izi] = port_inj_current(obj, x_, sysx, idx)
            % I = obj.port_inj_current(x_, sysx)
            % I = obj.port_inj_current(x_, sysx, idx)
            % [I, Iv1, Iv2] = obj.port_inj_current(...)
            % [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(...)
            % sysx : 1 = system x_, 0 = class aggregate x_

            if nargin < 4
                idx = [];
                if nargin < 3
                    sysx = 1;
                end
            end

            [Y, L, M, N, i, s] = obj.get_params(idx);
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear current injections and power injections
            nc = size(x_, 2);
            if nc > 1   %% Octave <= 5.2 has issues with broadcasting properly
                i = i * ones(1, nc);
                s = s * ones(1, nc);
            end
            if isempty(z_)
                Slin = M*v_ + s;
                I = Y*v_ + i + conj(Slin ./ vi_);
            else
                Slin = M*v_ + N*z_ + s;
                I = Y*v_ + L*z_ + i + conj(Slin ./ vi_);
            end

            if nargout > 1
                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                invdiagvic = sparse(1:ni, 1:ni, 1./conj(vi_), ni, ni);
                if any(any(M)) || any(any(N)) || any(s) || any(any(Y))
                    %% intermediate terms
                    if isempty(idx)     %% all ports
                        diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                    else                %% selected ports
                        diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                    end

                    [Iv1, Iv2] = obj.port_inj_current_jac( ...
                            n, v_, Y, M, invdiagvic, diagSlincJ);
                else
                    Iv1 = sparse(ni, n);
                    Iv2 = Iv1;
                end

                if nargout > 3
                    %% linear current term
                    Izr = L;
                    Izi = 1j * L;

                    %% + current from linear power term
                    IS_zr = invdiagvic * conj(N);   %% C
                    Izr = Izr + IS_zr;              %% +C
                    Izi = Izi - 1j * IS_zr;         %% -jC
                end

                if sysx
                    Ct = obj.C';
                    Iv1 = Iv1 * Ct;
                    Iv2 = Iv2 * Ct;
                    if nargout > 3
                        Dt = obj.D';
                        Izr = Izr * Dt;
                        Izi = Izi * Dt;
                    end
                end
            end

            %% general nonlinear current
            if ~isempty(obj.inln)
                Inln = cell(1, nargout);
                [Inln{:}] = obj.inln(x_, sysx, idx);
                I = I + Inln{1};
                if nargout > 1
                    Iv1 = Iv1 + Inln{2};
                    Iv2 = Iv2 + Inln{3};
                    if nargout > 3
                        Izr = Izr + Inln{4};
                        Izi = Izi + Inln{5};
                    end
                end
            elseif ~isempty(obj.snln)
                error('mp.form_ac/port_inj_current: Nonlinear current function not defined for corresponding nonlinear power function.')
            end
        end

        function [S, Sv1, Sv2, Szr, Szi] = port_inj_power(obj, x_, sysx, idx)
            % S = obj.port_inj_power(x_, sysx)
            % S = obj.port_inj_power(x_, sysx, idx)
            %
            % [S, Sv1, Sv2] = obj.port_inj_power(...)
            % [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(...)
            % sysx : 1 = system x_, 0 = class aggregate x_

            if nargin < 4
                idx = [];
                if nargin < 3
                    sysx = 1;
                end
            end

            [Y, L, M, N, i, s] = obj.get_params(idx);
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear current injections and power injections
            nc = size(x_, 2);
            if nc > 1   %% Octave <= 5.2 has issues with broadcasting properly
                i = i * ones(1, nc);
                s = s * ones(1, nc);
            end
            if isempty(z_)
                Ilin = Y*v_ + i;
                S = vi_ .* conj(Ilin) + M*v_ + s;
            else
                Ilin = Y*v_ + L*z_ + i;
                S = vi_ .* conj(Ilin) + M*v_ + N*z_ + s;
            end

            if nargout > 1
                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if any(any(Y)) || any(any(L)) || any(i) || any(any(M))
                    %% intermediate terms
                    diagv  = sparse(1:n, 1:n, v_, n, n);
                    if isempty(idx)     %% all ports
                        diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                    else                %% selected ports
                        diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                    end

                    [Sv1, Sv2] = obj.port_inj_power_jac( ...
                            n, v_, Y, M, diagv, diagvi, diagIlincJ);
                else                %% not a function of voltages
                    Sv1 = sparse(ni, n);
                    Sv2 = Sv1;
                end

                if nargout > 3
                    %% linear power term
                    Szr = N;
                    Szi = 1j * N;

                    %% + power from linear current term
                    SI_zr = diagvi * conj(L);   %% E
                    Szr = Szr + SI_zr;          %% +E
                    Szi = Szi - 1j * SI_zr;     %% -jE
                end

                if sysx
                    Ct = obj.C';
                    Sv1 = Sv1 * Ct;
                    Sv2 = Sv2 * Ct;
                    if nargout > 3
                        Dt = obj.D';
                        Szr = Szr * Dt;
                        Szi = Szi * Dt;
                    end
                end
            end

            %% general nonlinear power
            if ~isempty(obj.snln)
                Snln = cell(1, nargout);
                [Snln{:}] = obj.snln(x_, sysx, idx);
                S = S + Snln{1};
                if nargout > 1
                    Sv1 = Sv1 + Snln{2};
                    Sv2 = Sv2 + Snln{3};
                    if nargout > 3
                        Szr = Szr + Snln{4};
                        Szi = Szi + Snln{5};
                    end
                end
            elseif ~isempty(obj.inln)
                error('mp.form_ac/port_inj_power: Nonlinear power function not defined for corresponding nonlinear current function.')
            end
        end

        function H = port_inj_current_hess(obj, x_, lam, sysx, idx)
            % H = obj.port_inj_current_hess(x_, lam)
            % H = obj.port_inj_current_hess(x_, lam, sysx)
            % H = obj.port_inj_current_hess(x_, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear power injections
            if isempty(z_)
                Slin = M*v_ + s;
            else
                Slin = M*v_ + N*z_ + s;
            end

            nz = length(z_);
            n  = length(v_);    %% number of all port voltages
            if any(any(M)) || any(any(N)) || any(s) || any(any(Y))
                %% intermediate terms
                ni = length(vi_);   %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi_), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end

                [Iv1v1, Iv1v2, Iv2v2] = ...
                    obj.port_inj_current_hess_v(x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ);
                [Iv1zr, Iv1zi, Iv2zr, Iv2zi] = ...
                    obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, N, dlamJ);

                H = [   Iv1v1   Iv1v2   Iv1zr   Iv1zi;
                        Iv1v2.' Iv2v2   Iv2zr   Iv2zi;
                        Iv1zr.' Iv2zr.' sparse(nz, 2*nz);
                        Iv1zi.' Iv2zi.' sparse(nz, 2*nz)  ];
            else                %% not a function of voltages
                H = sparse(2*(n+nz), 2*(n+nz));
            end

            %% convert for system x_, if necessary
            if sysx
                C = obj.C;
                D = obj.D;
                [mc, nc] = size(C);
                [md, nd] = size(D);
                Ap = [  C sparse(mc, nc+2*nd);
                        sparse(mc,nc) C sparse(mc, 2*nd);
                        sparse(md, 2*nc) D sparse(md,nd);
                        sparse(md, 2*nc+nd) D ];
                H = Ap * H * Ap.';
            end

            %% general nonlinear current
            if ~isempty(obj.inln_hess)
                H = H + obj.inln_hess(x_, lam, sysx, idx);
            elseif ~isempty(obj.snln_hess)
                error('mp.form_ac/port_inj_current_hess: Nonlinear current Hessian not defined for corresponding nonlinear power Hessian.')
            end
        end

        function H = port_inj_power_hess(obj, x_, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x_, lam)
            % H = obj.port_inj_power_hess(x_, lam, sysx)
            % H = obj.port_inj_power_hess(x_, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear current injections
            if isempty(z_)
                Ilin = Y*v_ + i;
            else
                Ilin = Y*v_ + L*z_ + i;
            end

            nz = length(z_);
            n  = length(v_);    %% number of all port voltages
            if any(any(Y)) || any(any(L)) || any(i) || any(any(M))
                %% intermediate terms
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if isempty(idx)     %% all ports
                    diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end

                [Sv1v1, Sv1v2, Sv2v2] = ...
                    obj.port_inj_power_hess_v(x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ);
                [Sv1zr, Sv1zi, Sv2zr, Sv2zi] = ...
                    obj.port_inj_power_hess_vz(x_, lam, v_, z_, diagvi, L, dlamJ);

                H = [   Sv1v1   Sv1v2   Sv1zr   Sv1zi;
                        Sv1v2.' Sv2v2   Sv2zr   Sv2zi;
                        Sv1zr.' Sv2zr.' sparse(nz, 2*nz);
                        Sv1zi.' Sv2zi.' sparse(nz, 2*nz)  ];
            else                %% not a function of voltages
                H = sparse(2*(n+nz), 2*(n+nz));
            end

            %% convert for system x_, if necessary
            if sysx
                C = obj.C;
                D = obj.D;
                [mc, nc] = size(C);
                [md, nd] = size(D);
                Ap = [  C sparse(mc, nc+2*nd);
                        sparse(mc,nc) C sparse(mc, 2*nd);
                        sparse(md, 2*nc) D sparse(md,nd);
                        sparse(md, 2*nc+nd) D ];
                H = Ap * H * Ap.';
            end

            %% general nonlinear power
            if ~isempty(obj.snln_hess)
                H = H + obj.snln_hess(x_, lam, sysx, idx);
            elseif ~isempty(obj.inln_hess)
                error('mp.form_ac/port_inj_power_hess: Nonlinear power Hessian not defined for corresponding nonlinear current Hessian.')
            end
        end

        function [h, dh] = port_apparent_power_lim_fcn(obj, x_, nm, idx, hmax)
            %% branch squared apparent power flow constraints

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                Sx = [Sv1 Sv2 Szr Szi];
                n = length(S);
                dS = spdiags(S, 0, n, n);
                dh = 2 * (real(dS) * real(Sx) + imag(dS) * imag(Sx));
            else
                S = obj.port_inj_power(x_, 1, idx);
            end
            h = conj(S) .* S - hmax;
        end

        function [h, dh] = port_active_power_lim_fcn(obj, x_, nm, idx, hmax)
            %% branch active power flow constraints

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                P = real(S);
                dh = real([Sv1 Sv2 Szr Szi]);
            else
                P = real(obj.port_inj_power(x_, 1, idx));
            end
            h = P - hmax;
        end

        function [h, dh] = port_active_power2_lim_fcn(obj, x_, nm, idx, hmax)
            %% branch squared active power flow constraints

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                Sx = [Sv1 Sv2 Szr Szi];
                P = real(S);
                n = length(S);
                dP = spdiags(P, 0, n, n);
                dh = 2 * real(dP) * real(Sx);
            else
                P = real(obj.port_inj_power(x_, 1, idx));
            end
            h = P .* P - hmax;
        end

        function [h, dh] = port_current_lim_fcn(obj, x_, nm, idx, hmax)
            %% branch squared current constraints

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1, idx);
                Ix = [Iv1 Iv2 Izr Izi];
                n = length(I);
                dI = spdiags(I, 0, n, n);
                dh = 2 * (real(dI) * real(Ix) + imag(dI) * imag(Ix));
            else
                I = obj.port_inj_power(x_, 1, idx);
            end
            h = conj(I) .* I - hmax;
        end

        function d2H = port_apparent_power_lim_hess(obj, x_, lam, nm, idx)
            %% branch squared apparent power flow Hessian
            nlam = length(lam);

            [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
            n = length(S);
            dSc = spdiags(conj(S), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Sx = [Sv1 Sv2 Szr Szi];

            d2H = 2 * real( obj.port_inj_power_hess(x_, dSc * lam, 1, idx) + ...
                        Sx.' * dlam * conj(Sx) );
        end

        function d2H = port_active_power_lim_hess(obj, x_, lam, nm, idx)
            %% branch active power flow Hessian

            d2H = real( obj.port_inj_power_hess(x_, lam, 1, idx) );
        end

        function d2H = port_active_power2_lim_hess(obj, x_, lam, nm, idx)
            %% branch squared active power flow Hessian
            nlam = length(lam);

            [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
            n = length(S);
            dP = spdiags(real(S), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Px = real([Sv1 Sv2 Szr Szi]);

            d2H = 2 * real( obj.port_inj_power_hess(x_, dP * lam, 1, idx) + ...
                        Px.' * dlam * Px );
        end

        function d2H = port_current_lim_hess(obj, x_, lam, nm, idx)
            %% branch squared current Hessian
            nlam = length(lam);

            [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1, idx);
            n = length(I);
            dIc = spdiags(conj(I), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Ix = [Iv1 Iv2 Izr Izi];

            d2H = 2 * real( obj.port_inj_power_hess(x_, dIc * lam, 1, idx) + ...
                        Ix.' * dlam * conj(Ix) );
        end
    end     %% methods
end         %% classdef
