classdef opf_model < opt_model
%OPF_MODEL  Constructor for OPF model class.
%   OM = OPF_MODEL(MPC)
%
%   This class implements the OPF model object used to encapsulate
%   a given OPF problem formulation. It allows for access to optimization
%   variables, constraints and costs in named blocks, keeping track of the
%   ordering and indexing of the blocks as variables, constraints and costs
%   are added to the problem.
%
%   This class is a sub-class of OPT_MODEL that adds the 'mpc'
%   field for storing the MATPOWER case struct used to build the object
%   along with the get_mpc() method.
%
%   It also add the 'cost' field and the following three methods for
%   implementing the legacy user-defined OPF costs:
%       add_legacy_cost
%       params_legacy_cost
%       eval_legacy_cost
%
%   The following deprecated OPT_MODEL methods, from an older OPF_MODEL
%   implementation have also been moved back to OPF_MODEL:
%       add_costs (use add_quad_cost, add_nln_cost or add_legacy_cost)
%       add_constraints (use add_lin_constraint or add_nln_constraint)
%       add_vars (use add_var)
%       getv (use params_var)
%       linear_constraints (use params_lin_constraint)
%       build_cost_params (use params_legacy_cost)
%       compute_cost (use eval_legacy_cost)
%
%   The following is the structure of the data in the OPF model object.
%
%   om
%       <opt_model fields> - see OPT_MODEL for details
%       .cost       - data for legacy user-defined costs
%           .idx
%               .i1 - starting row index within full N matrix
%               .iN - ending row index within full N matrix
%               .N  - number of rows in this cost block in full N matrix
%           .N      - total number of rows in full N matrix
%           .NS     - number of cost blocks
%           .data   - data for each user-defined cost block
%               .N  - see help for ADD_LEGACY_COST for details
%               .H  -               "
%               .Cw -               "
%               .dd -               "
%               .rr -               "
%               .kk -               "
%               .mm -               "
%               .vs - cell array of variable sets that define xx for this
%                     cost block, where the N for this block multiplies xx
%           .order  - struct array of names/indices for cost blocks in the
%                     order they appear in the rows of the full N matrix
%               .name   - name of the block, e.g. R
%               .idx    - indices for name, {2,3} => R(2,3)
%       .mpc        - MATPOWER case struct used to create this model object
%           .baseMVA
%           .bus
%           .branch
%           .gen
%           .gencost
%           .A  (if present, must have l, u)
%           .l
%           .u
%           .N  (if present, must have fparm, H, Cw)
%           .fparm
%           .H
%           .Cw
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost = struct( ...
            'idx', struct( ...
                'i1', struct(), ...
                'iN', struct(), ...
                'N', struct() ), ...
            'N', 0, ...
            'NS', 0, ...
            'order', struct( ...
                'name', [], ...
                'idx', [] ), ...
            'data', struct( ...
                'N', struct(), ...
                'H', struct(), ...
                'Cw', struct(), ...
                'dd', struct(), ...
                'rh', struct(), ...
                'kk', struct(), ...
                'mm', struct(), ...
                'vs', struct() ), ...
            'params', [] );
        mpc = struct();
    end

    methods
        %% constructor
        function om = opf_model(mpc)
            args = {};
            have_mpc = 0;
            if nargin > 0
                if isa(mpc, 'opf_model')
                    args = { mpc };
                elseif isstruct(mpc)
                    have_mpc = 1;
                end
            end

            om@opt_model(args{:});
            
            if have_mpc
                om.mpc = mpc;
            end
        end
    end
end
