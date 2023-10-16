classdef form_dc < mp.form
% mp.form_dc - Base class for |MATPOWER| DC **formulations**.
%
% Used as a mix-in class for all **network model element** classes
% with a DC network model formulation. That is, each concrete network model
% element class with a DC formulation must inherit, at least indirectly, from
% both mp.nm_element and mp.form_dc.
%
% mp.form_dc defines the port active power injection as a linear function
% of the state variables :math:`\x`, that is, the voltage angles :math:`\Va`
% and non-voltage states :math:`\z`, as
%
% .. math::
%
%    \gP(\x) &= \left[\begin{array}{cc}\BB & \KK\end{array}\right] \x + \pv \\
%    &= \BB \Va + \KK \z + \pv,
%
% where :math:`\BB`, :math:`\KK`, and :math:`\pv` are the model parameters.
%
% For more details, see the :ref:`sec_nm_formulations_dc` section in the
% |MATPOWER-Dev-Manual| and the derivations in |TN5|.
%
% mp.form_dc Properties:
%   * B - :math:`n_p n_k \times n_n` matrix :math:`\BB` of model parameters
%   * K - :math:`n_p n_k \times n_z` matrix :math:`\KK` of model parameters
%   * p - :math:`n_p n_k \times 1` vector :math:`\pv` of model parameters
%   * params_ncols - specify number of columns for each parameter
%
% mp.form_dc Methods:
%   * form_name - get char array w/name of formulation (``'DC'``)
%   * form_tag - get char array w/short label of formulation (``'dc'``)
%   * model_params - get network model element parameters (``{'B', 'K', 'p'}``)
%   * model_vvars - get cell array of names of voltage state variables (``{'va'}``)
%   * model_zvars - get cell array of names of non-voltage state variables (``{'z'}``)
%   * port_inj_power - compute port power injections from network state
%
% See also mp.form, mp.form_ac, mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % *(double)* :math:`n_p n_k \times n_n` matrix :math:`\BB`
        % of model parameter coefficients for :math:`\Va`
        B = [];

        % *(double)* :math:`n_p n_k \times n_z` matrix :math:`\KK`
        % of model parameter coefficients for :math:`\z`
        K = [];
        p = []; % *(double)* :math:`n_p n_k \times 1` vector :math:`\pv` of model parameters

        % *(struct)* specify number of columns for each parameter,
        % where
        %
        %   - 1 => single column (i.e. a vector)
        %   - 2 => :math:`n_p` columns
        %   - 3 => :math:`n_z` columns
        param_ncols = struct('B', 2, 'K', 3, 'p', 1);   
    end

    methods
        function name = form_name(obj)
            % Get user-readable name of formulation, i.e. ``'DC'``.
            %
            % See :meth:`mp.form.form_name`.

            name = 'DC';
        end

        function tag = form_tag(obj)
            % Get short label of formulation, i.e. ``'dc'``.
            %
            % See :meth:`mp.form.form_tag`.

            tag = 'dc';
        end

        function params = model_params(obj)
            % Get cell array of names of model parameters, i.e. ``{'B', 'K', 'p'}``.
            %
            % See :meth:`mp.form.model_params`.

           params = {'B', 'K', 'p'};
        end

        function vtypes = model_vvars(obj)
            % Get cell array of names of voltage state variables, i.e. ``{'va'}``.
            %
            % See :meth:`mp.form.model_vvars`.

            vtypes = {'va'};
        end

        function vtypes = model_zvars(obj)
            % Get cell array of names of non-voltage state variables, i.e. ``{'z'}``.
            %
            % See :meth:`mp.form.model_zvars`.

            vtypes = {'z'};
        end

        function P = port_inj_power(obj, x, sysx, idx)
            % Compute port power injections from network state.
            % ::
            %
            %   P = nme.port_inj_power(x, sysx, idx)
            %
            % Compute the active power injections for all or a selected
            % subset of ports.
            %
            % The state can be provided as a stacked aggregate of the state
            % variables (port voltages and non-voltage states) for the full
            % collection of network model elements of this type, or as the
            % combined state for the entire network.
            %
            % Inputs:
            %   x (double) : state vector :math:`\x`
            %   sysx (0 or 1) : which state is provided in ``x``
            %
            %       - 0 -- class aggregate state
            %       - 1 -- *(default)* full system state
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns injections
            %       corresponding to all ports
            %
            % Outputs:
            %   P (double) : vector of port power injections, :math:`\gP(\x)`

            if nargin < 4
                idx = [];
                if nargin < 3
                    sysx = 1;
                end
            end

            [B, K, p] = obj.get_params(idx);
            [v, z] = obj.x2vz(x, sysx);

            if isempty(z)
                P = B*v + p;
            else
                P = B*v + K*z + p;
            end
%             if sysx
%                 Ct = obj.C';
%                 Dt = obj.D';
%                 if nargin < 4       %% all ports
%                     P = B*Ct*v + K*Dt*z + p;
%                 else                %% selected ports
%                     P = B(idx, :)*Ct*v + K(idx, :)*Dt*z + p(idx);
%                 end
%             else
%                 if nargin < 4       %% all ports
%                     P = B*v + K*z + p;
%                 else                %% selected ports
%                     P = B(idx, :)*v + K(idx, :)*z + p(idx);
%                 end
%             end
        end
    end     %% methods
end         %% classdef
