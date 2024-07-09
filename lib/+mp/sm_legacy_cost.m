classdef sm_legacy_cost < mp.set_manager_opt_model
% mp.sm_legacy_cost -  MP Set Manager class for legacy costs.
% ::
%
%   sm = mp.sm_legacy_cost()
%   sm = mp.sm_legacy_cost(label)
%
% MP Set Manager class for legacy costs of the form described in
% opf_model.add_legacy_cost.
%
% mp.sm_legacy_cost Properties:
%   * cache - struct for caching aggregated parameters for the set
%
% mp.sm_legacy_cost Methods:
%   * sm_legacy_cost - constructor
%
% See also mp.set_manager, mp.set_manager_opt_model.

%   MATPOWER
%   Copyright (c) 2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % struct for caching aggregated parameters for legacy costs
        cache = [];
    end     %% properties

    methods
        function obj = sm_legacy_cost(varargin)
            % Constructor.
            % ::
            %
            %   sm = mp.sm_legacy_cost(label)

            es = struct();  %% empty struct
            obj@mp.set_manager_opt_model(varargin{:});
            obj.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
        end
    end     %% methods
end         %% classdef
