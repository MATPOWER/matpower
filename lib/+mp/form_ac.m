classdef (Abstract) form_ac < mp.form
% mp.form_ac - Abstract base class for |MATPOWER| AC **formulations**.
%
% Used as a mix-in class for all **network model element** classes
% with an AC network model formulation. That is, each concrete network model
% element class with an AC formulation must inherit, at least indirectly, from
% both mp.nm_element and mp.form_ac.
%
% mp.form_ac defines the complex port injections as functions of the state
% variables :math:`\X`, that is, the complex voltages :math:`\V` and
% non-voltage states :math:`\Z`. They are defined in terms of 3 components,
% the linear current injection and linear power injection components,
% 
% .. math::
% 
%    \Ilin(\X) &= \left[\begin{array}{cc}\YY & \LL\end{array}\right] \X + \iv \\
%    &= \YY \V + \LL \Z + \iv
% 
% .. math::
% 
%    \Slin(\X) &= \left[\begin{array}{cc}\MM & \NN\end{array}\right] \X + \sv \\
%    &= \MM \V + \NN \Z + \sv,
%
% and an arbitrary nonlinear injection component represented by
% :math:`\Snln(\X)` or :math:`\Inln(\X)`. The full complex power and current
% port injection functions implemented by mp.form_ac, are respectively 
%
% .. math::
%
%     \GS(\X) &= \dV \conj{\left( \Ilin(\X) \right)} + \Slin(\X) + \Snln(\X) \\
%     &= \dV \conj{\left( \YY \V + \LL \Z + \iv \right)} + \MM \V + \NN \Z + \sv + \Snln(\X)
%
% .. math::
%
%   \GI(\X) &= \Ilin(\X) + \cdiag{\Slin(\X)} \inVc + \Inln(\X) \\
%   &= \YY \V + \LL \Z + \iv + \cdiag{\MM \V + \NN \Z + \sv} \inVc + \Inln(\X)
%
% where :math:`\YY`, :math:`\LL`, :math:`\MM`, :math:`\NN`, :math:`\iv`, and
% :math:`\sv`, along with :math:`\Snln(\X)` or :math:`\Inln(\X)`, are the
% model parameters.
%
% For more details, see the :ref:`sec_nm_formulations_ac` section in the
% |MATPOWER-Dev-Manual| and the derivations in |TN5|.
%
% mp.form_dc Properties:
%   * Y - :math:`n_p n_k \times n_n` matrix :math:`\YY` of model parameters
%   * L - :math:`n_p n_k \times n_\Z` matrix :math:`\LL` of model parameters
%   * M - :math:`n_p n_k \times n_n` matrix :math:`\MM` of model parameters
%   * N - :math:`n_p n_k \times n_\Z` matrix :math:`\NN` of model parameters
%   * i - :math:`n_p n_k \times 1` vector :math:`\iv` of model parameters
%   * s - :math:`n_p n_k \times 1` vector :math:`\sv` of model parameters
%   * params_ncols - specify number of columns for each parameter
%   * inln - function to compute :math:`\Inln(\X)`
%   * snln - function to compute :math:`\Snln(\X)`
%   * inln_hess - function to compute Hessian of :math:`\Inln(\X)`
%   * snln_hess - function to compute Hessian of :math:`\Snln(\X)`
%
% mp.form_dc Methods:
%   * model_params - get network model element parameters (``{'Y', 'L', 'M', 'N', 'i', 's'}``)
%   * model_zvars - get cell array of names of non-voltage state variables (``{'zr', 'zi'}``)
%   * port_inj_current - compute port current injections from network state
%   * port_inj_power - compute port power injections from network state
%   * port_inj_current_hess - compute Hessian of port current injections
%   * port_inj_power_hess - compute Hessian of port power injections
%   * port_inj_current_jac - abstract method to compute voltage-related Jacobian terms
%   * port_inj_current_hess_v - abstract method to compute voltage-related Hessian terms
%   * port_inj_current_hess_vz - abstract method to compute voltage-related Hessian terms
%   * port_inj_power_jac - abstract method to compute voltage-related Jacobian terms
%   * port_inj_power_hess_v - abstract method to compute voltage-related Hessian terms
%   * port_inj_power_hess_vz - abstract method to compute voltage-related Hessian terms
%   * port_apparent_power_lim_fcn - compute port squared apparent power injection constraints
%   * port_active_power_lim_fcn - compute port active power injection constraints
%   * port_active_power2_lim_fcn - compute port squared active power injection constraints
%   * port_current_lim_fcn - compute port squared current injection constraints
%   * port_apparent_power_lim_hess - compute port squared apparent power injection Hessian
%   * port_active_power_lim_hess - compute port active power injection Hessian
%   * port_active_power2_lim_hess - compute port squared active power injection Hessian
%   * port_current_lim_hess - compute port squared current injection Hessian
%   * aux_data_va_vm - abstract method to return voltage angles/magnitudes from auxiliary data
%
% See also mp.form, mp.form_acc, mp.form_acp, mp.form_dc, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % *(double)* :math:`n_p n_k \times n_n` matrix :math:`\YY`
        % of model parameter coefficients for :math:`\V`
        Y = [];

        % *(double)* :math:`n_p n_k \times n_\Z` matrix :math:`\LL`
        % of model parameter coefficients for :math:`\Z`
        L = [];

        % *(double)* :math:`n_p n_k \times n_n` matrix :math:`\MM`
        % of model parameter coefficients for :math:`\V`
        M = [];

        % *(double)* :math:`n_p n_k \times n_\Z` matrix :math:`\NN`
        % of model parameter coefficients for :math:`\Z`
        N = [];
        i = []; % *(double)* :math:`n_p n_k \times 1` vector :math:`\iv` of model parameters
        s = []; % *(double)* :math:`n_p n_k \times 1` vector :math:`\sv` of model parameters

        % *(struct)* specify number of columns for each parameter,
        % where
        %
        %   - 1 => single column (i.e. a vector)
        %   - 2 => :math:`n_p` columns
        %   - 3 => :math:`n_\Z` columns
        param_ncols = struct('Y', 2, 'L', 3, 'M', 2, 'N', 3, 'i', 1, 's', 1);
        inln = '';      % *(function handle)* function to compute :math:`\Inln(\X)`
        snln = '';      % *(function handle)* function to compute :math:`\Snln(\X)`
        inln_hess = ''; % *(function handle)* function to compute Hessian of :math:`\Inln(\X)`
        snln_hess = ''; % *(function handle)* function to compute Hessian of :math:`\Snln(\X)`
    end

    methods
        function params = model_params(obj)
            % Get cell array of names of model parameters, i.e. ``{'Y', 'L', 'M', 'N', 'i', 's'}``.
            %
            % See :meth:`mp.form.model_params`.

           params = {'Y', 'L', 'M', 'N', 'i', 's'};
        end

        function vtypes = model_zvars(obj)
            % Get cell array of names of non-voltage state variables, i.e. ``{'zr', 'zi'}``.
            %
            % See :meth:`mp.form.model_zvars`.

            vtypes = {'zr', 'zi'};
        end

        function [I, Iv1, Iv2, Izr, Izi] = port_inj_current(obj, x_, sysx, idx)
            % Compute port complex current injections from network state.
            % ::
            %
            %   I = nme.port_inj_current(x_, sysx)
            %   I = nme.port_inj_current(x_, sysx, idx)
            %   [I, Iv1, Iv2] = nme.port_inj_current(...)
            %   [I, Iv1, Iv2, Izr, Izi] = nme.port_inj_current(...)
            %
            % Compute the complex current injections for all or a selected
            % subset of ports and, optionally, the components of the Jacobian,
            % that is, the sparse matrices of partial derivatives with respect
            % to each real component of the state. The voltage portion, which
            % depends on the formulation (polar vs cartesian), is delegated
            % to the  :meth:`port_inj_current_jac` method implemented by the
            % appropriate subclass.
            %
            % The state can be provided as a stacked aggregate of the state
            % variables (port voltages and non-voltage states) for the full
            % collection of network model elements of this type, or as the
            % combined state for the entire network.
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   sysx (0 or 1) : which state is provided in ``x_``
            %
            %       - 0 -- class aggregate state
            %       - 1 -- *(default)* full system state
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Outputs:
            %   I (complex double) : vector of port complex current injections,
            %       :math:`\GI(\X)`
            %   Iv1 (complex double) : Jacobian of port complex current
            %       injections w.r.t 1st voltage component,
            %       :math:`\GI_\Va` (polar) or :math:`\GI_\Vr`
            %       (cartesian)
            %   Iv2 (complex double) : Jacobian of port complex current
            %       injections w.r.t 2nd voltage component,
            %       :math:`\GI_\Vm` (polar) or :math:`\GI_\Vi`
            %       (cartesian)
            %   Izr (complex double) : Jacobian of port complex current
            %       injections w.r.t real part of non-voltage state,
            %       :math:`\GI_\Zr`
            %   Izi (complex double) : Jacobian of port complex current
            %       injections w.r.t imaginary part of non-voltage state,
            %       :math:`\GI_\Zi`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power.

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
                error('mp.form_ac.port_inj_current: Nonlinear current function not defined for corresponding nonlinear power function.')
            end
        end

        function [S, Sv1, Sv2, Szr, Szi] = port_inj_power(obj, x_, sysx, idx)
            % Compute port complex power injections from network state.
            % ::
            %
            %   S = nme.port_inj_power(x_, sysx)
            %   S = nme.port_inj_power(x_, sysx, idx)
            %   [S, Sv1, Sv2] = nme.port_inj_power(...)
            %   [S, Sv1, Sv2, Szr, Szi] = nme.port_inj_power(...)
            %
            % Compute the complex power injections for all or a selected
            % subset of ports and, optionally, the components of the Jacobian,
            % that is, the sparse matrices of partial derivatives with respect
            % to each real component of the state. The voltage portion, which
            % depends on the formulation (polar vs cartesian), is delegated
            % to the  :meth:`port_inj_power_jac` method implemented by the
            % appropriate subclass.
            %
            % The state can be provided as a stacked aggregate of the state
            % variables (port voltages and non-voltage states) for the full
            % collection of network model elements of this type, or as the
            % combined state for the entire network.
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   sysx (0 or 1) : which state is provided in ``x_``
            %
            %       - 0 -- class aggregate state
            %       - 1 -- *(default)* full system state
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Outputs:
            %   S (complex double) : vector of port complex power injections,
            %       :math:`\GS(\X)`
            %   Sv1 (complex double) : Jacobian of port complex power
            %       injections w.r.t 1st voltage component,
            %       :math:`\GS_\Va` (polar) or :math:`\GS_\Vr`
            %       (cartesian)
            %   Sv2 (complex double) : Jacobian of port complex power
            %       injections w.r.t 2nd voltage component,
            %       :math:`\GS_\Vm` (polar) or :math:`\GS_\Vi`
            %       (cartesian)
            %   Szr (complex double) : Jacobian of port complex power
            %       injections w.r.t real part of non-voltage state,
            %       :math:`\GS_\Zr`
            %   Szi (complex double) : Jacobian of port complex power
            %       injections w.r.t imaginary part of non-voltage state,
            %       :math:`\GS_\Zi`.
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_current.

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
                error('mp.form_ac.port_inj_power: Nonlinear power function not defined for corresponding nonlinear current function.')
            end
        end

        function H = port_inj_current_hess(obj, x_, lam, sysx, idx)
            % Compute Hessian of port current injections from network state.
            % ::
            %
            %   H = nme.port_inj_current_hess(x_, lam)
            %   H = nme.port_inj_current_hess(x_, lam, sysx)
            %   H = nme.port_inj_current_hess(x_, lam, sysx, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the port current injection Jacobian by a vector :math:`\lam`.
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\lam` of multipliers, one for each port
            %   sysx (0 or 1) : which state is provided in ``x_``
            %
            %       - 0 -- class aggregate state
            %       - 1 -- *(default)* full system state
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Outputs:
            %   H (complex double) : sparse Hessian matrix of port complex
            %       current injections corresponding to specified
            %       :math:`\lam`, namely :math:`\GI_{\X\X}(\lam)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_current.

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
                error('mp.form_ac.port_inj_current_hess: Nonlinear current Hessian not defined for corresponding nonlinear power Hessian.')
            end
        end

        function H = port_inj_power_hess(obj, x_, lam, sysx, idx)
            % Compute Hessian of port power injections from network state.
            % ::
            %
            %   H = nme.port_inj_power_hess(x_, lam)
            %   H = nme.port_inj_power_hess(x_, lam, sysx)
            %   H = nme.port_inj_power_hess(x_, lam, sysx, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the port power injection Jacobian by a vector :math:`\lam`.
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\lam` of multipliers, one for each port
            %   sysx (0 or 1) : which state is provided in ``x_``
            %
            %       - 0 -- class aggregate state
            %       - 1 -- *(default)* full system state
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Outputs:
            %   H (complex double) : sparse Hessian matrix of port complex
            %       power injections corresponding to specified
            %       :math:`\lam`, namely :math:`\GS_{\X\X}(\lam)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power.

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
                error('mp.form_ac.port_inj_power_hess: Nonlinear power Hessian not defined for corresponding nonlinear current Hessian.')
            end
        end

        function [Iv1, Iv2] = port_inj_current_jac(obj, n, v_, Y, M, invdiagvic, diagSlincJ)
            % Abstract method to compute voltage-related Jacobian terms.
            %
            % Called by port_inj_current() to compute voltage-related Jacobian
            % terms. See mp.form_acc.port_inj_current_jac and
            % mp.form_acp.port_inj_current_jac for details.

            error('mp.form_ac.port_inj_current_jac: must be implemented in subclass');
        end

        function [Iv1v1, Iv1v2, Iv2v2] = port_inj_current_hess_v(obj, x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            % Abstract method to compute voltage-related Hessian terms.
            %
            % Called by port_inj_current_hess() to compute voltage-related Hessian
            % terms. See mp.form_acc.port_inj_current_hess_v and
            % mp.form_acp.port_inj_current_hess_v for details.

            error('mp.form_ac.port_inj_current_hess_v: must be implemented in subclass');
        end

        function [Iv1zr, Iv1zi, Iv2zr, Iv2zi] = port_inj_current_hess_vz(obj, x_, lam, v_, z_, diaginvic, N, dlamJ)
            % Abstract method to compute voltage-related Hessian terms.
            %
            % Called by port_inj_current_hess() to compute
            % voltage/non-voltage-related Hessian terms. See
            % mp.form_acc.port_inj_current_hess_vz and
            % mp.form_acp.port_inj_current_hess_vz for details.

            error('mp.form_ac.port_inj_current_hess_vz: must be implemented in subclass');
        end

        function [Sv1, Sv2] = port_inj_power_jac(obj, n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % Abstract method to compute voltage-related Jacobian terms.
            %
            % Called by port_inj_power() to compute voltage-related Jacobian
            % terms. See mp.form_acc.port_inj_power_jac and
            % mp.form_acp.port_inj_power_jac for details.

            error('mp.form_ac.port_inj_power_jac: must be implemented in subclass');
        end

        function [Sv1v1, Sv1v2, Sv2v2] = port_inj_power_hess_v(obj, x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            % Abstract method to compute voltage-related Hessian terms.
            %
            % Called by port_inj_power_hess() to compute voltage-related Hessian
            % terms. See mp.form_acc.port_inj_power_hess_v and
            % mp.form_acp.port_inj_power_hess_v for details.

            error('mp.form_ac.port_inj_power_hess_v: must be implemented in subclass');
        end

        function [Sv1zr, Sv1zi, Sv2zr, Sv2zi] = port_inj_power_hess_vz(obj, x_, lam, v_, z_, diagvi, L, dlamJ)
            % Abstract method to compute voltage-related Hessian terms.
            %
            % Called by port_inj_power_hess() to compute
            % voltage/non-voltage-related Hessian terms. See
            % mp.form_acc.port_inj_power_hess_vz and
            % mp.form_acp.port_inj_power_hess_vz for details.

            error('mp.form_ac.port_inj_power_hess_vz: must be implemented in subclass');
        end

        function [h, dh] = port_apparent_power_lim_fcn(obj, x_, nm, idx, hmax)
            % Compute port squared apparent power injection constraints.
            % ::
            %
            %   h = nme.port_apparent_power_lim_fcn(x_, nm, idx, hmax)
            %   [h, dh] = nme.port_apparent_power_lim_fcn(x_, nm, idx, hmax)
            %
            % Compute constraint function and optionally the Jacobian for
            % the limit on port squared apparent power injections based on
            % complex outputs of port_inj_power().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %   hmax (double) : vector of squared apparent power limits
            %
            % Outputs:
            %   h (double) : constraint function, :math:`\h^\mathrm{flow}(\x)`
            %   dh (double) : constraint Jacobian, :math:`\h_\x^\mathrm{flow}`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power.

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
            % Compute port active power injection constraints.
            % ::
            %
            %   h = nme.port_active_power_lim_fcn(x_, nm, idx, hmax)
            %   [h, dh] = nme.port_active_power_lim_fcn(x_, nm, idx, hmax)
            %
            % Compute constraint function and optionally the Jacobian for
            % the limit on port active power injections based on
            % complex outputs of port_inj_power().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %   hmax (double) : vector of active power limits
            %
            % Outputs:
            %   h (double) : constraint function, :math:`\h^\mathrm{flow}(\x)`
            %   dh (double) : constraint Jacobian, :math:`\h_\x^\mathrm{flow}`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power.

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
            % Compute port squared active power injection constraints.
            % ::
            %
            %   h = nme.port_active_power2_lim_fcn(x_, nm, idx, hmax)
            %   [h, dh] = nme.port_active_power2_lim_fcn(x_, nm, idx, hmax)
            %
            % Compute constraint function and optionally the Jacobian for
            % the limit on port squared active power injections based on
            % complex outputs of port_inj_power().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %   hmax (double) : vector of squared active power limits
            %
            % Outputs:
            %   h (double) : constraint function, :math:`\h^\mathrm{flow}(\x)`
            %   dh (double) : constraint Jacobian, :math:`\h_\x^\mathrm{flow}`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power.

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
            % Compute port squared current injection constraints.
            % ::
            %
            %   h = nme.port_current_lim_fcn(x_, nm, idx, hmax)
            %   [h, dh] = nme.port_current_lim_fcn(x_, nm, idx, hmax)
            %
            % Compute constraint function and optionally the Jacobian for
            % the limit on port squared current injections based on
            % complex outputs of port_inj_current().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %   hmax (double) : vector of squared current limits
            %
            % Outputs:
            %   h (double) : constraint function, :math:`\h^\mathrm{flow}(\x)`
            %   dh (double) : constraint Jacobian, :math:`\h_\x^\mathrm{flow}`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_current.

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
            % Compute port squared apparent power injection Hessian.
            % ::
            %
            %   d2H = nme.port_apparent_power_lim_hess(x_, lam, nm, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the constraint Jacobian by a vector :math:`\muv`. Results are
            % based on the complex outputs of port_inj_power() and
            % port_inj_power_hess().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\muv` of multipliers, one for each port
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Output:
            %   d2H (double) : sparse constraint Hessian matrix, :math:`\h_{\x\x}^\mathrm{flow}(\muv)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power, port_inj_power_hess.

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
            % Compute port active power injection Hessian.
            % ::
            %
            %   d2H = nme.port_active_power_lim_hess(x_, lam, nm, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the constraint Jacobian by a vector :math:`\muv`. Results are
            % based on the complex outputs of port_inj_power() and
            % port_inj_power_hess().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\muv` of multipliers, one for each port
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Output:
            %   d2H (double) : sparse constraint Hessian matrix, :math:`\h_{\x\x}^\mathrm{flow}(\muv)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power, port_inj_power_hess.

            d2H = real( obj.port_inj_power_hess(x_, lam, 1, idx) );
        end

        function d2H = port_active_power2_lim_hess(obj, x_, lam, nm, idx)
            % Compute port squared active power injection Hessian.
            % ::
            %
            %   d2H = nme.port_active_power2_lim_hess(x_, lam, nm, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the constraint Jacobian by a vector :math:`\muv`. Results are
            % based on the complex outputs of port_inj_power() and
            % port_inj_power_hess().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\muv` of multipliers, one for each port
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Output:
            %   d2H (double) : sparse constraint Hessian matrix, :math:`\h_{\x\x}^\mathrm{flow}(\muv)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_power, port_inj_power_hess.

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
            % Compute port squared current injection Hessian.
            % ::
            %
            %   d2H = nme.port_current_lim_hess(x_, lam, nm, idx)
            %
            % Compute a sparse Hessian matrix for all or a selected subset of
            % ports. Rather than a full, 3-dimensional Hessian, it computes the
            % Jacobian of the vector obtained by muliplying the transpose of
            % the constraint Jacobian by a vector :math:`\muv`. Results are
            % based on the complex outputs of port_inj_current() and
            % port_inj_current_hess().
            %
            % Inputs:
            %   x_ (complex double) : state vector :math:`\X`
            %   lam (double) : vector :math:`\muv` of multipliers, one for each port
            %   nm (mp.net_model) : network model object
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns results
            %       corresponding to all ports
            %
            % Output:
            %   d2H (double) : sparse constraint Hessian matrix, :math:`\h_{\x\x}^\mathrm{flow}(\muv)`
            %
            % For details on the derivations of the formulas used, see |TN5|.
            %
            % See also port_inj_current, port_inj_current_hess.

            nlam = length(lam);

            [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1, idx);
            n = length(I);
            dIc = spdiags(conj(I), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Ix = [Iv1 Iv2 Izr Izi];

            d2H = 2 * real( obj.port_inj_power_hess(x_, dIc * lam, 1, idx) + ...
                        Ix.' * dlam * conj(Ix) );
        end

        function [va, vm] = aux_data_va_vm(obj, ad)
            % Abstract method to return voltage angles/magnitudes from auxiliary data.
            % ::
            %
            %   [va, vm] = nme.aux_data_va_vm(ad)
            %
            % Input:
            %   ad (struct) : struct of auxiliary data
            %
            % Outputs:
            %   va (double) : vector of voltage angles corresponding to
            %       voltage information stored in auxiliary data
            %   vm (double) : vector of voltage magnitudes corresponding to
            %       voltage information stored in auxiliary data
            %
            % Implemented by mp.form_acc.aux_data_va_vm and
            % mp.form_acp.aux_data_va_vm.

            error('mp.form_ac.aux_data_va_vm: must be implemented in subclass');
        end
    end     %% methods
end         %% classdef
