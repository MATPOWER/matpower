classdef opf_model < opt_model
% opf_model - Legacy |MATPOWER| OPF model class.
% ::
%
%   OM = OPF_MODEL(MPC)
%
%   This class implements the OPF model object used to encapsulate
%   a given OPF problem formulation. It allows for access to optimization
%   variables, constraints and costs in named blocks, keeping track of the
%   ordering and indexing of the blocks as variables, constraints and costs
%   are added to the problem.
%
%   This class is a subclass of OPT_MODEL that adds the 'mpc'
%   field for storing the MATPOWER case struct used to build the object
%   along with the get_mpc() method.
%
%   It also adds the 'cost' field and the following three methods for
%   implementing the legacy user-defined OPF costs:
%       add_legacy_cost
%       params_legacy_cost
%       eval_legacy_cost
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
% See also opt_model.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost = [];          % data for legacy user-defined costs
        mpc = struct();     % |MATPOWER| case struct from which ``om`` was built
    end     %% properties

    methods
        function om = opf_model(mpc)
            % Constructor.
            % ::
            %
            %   om = opf_model()
            %   om = opf_model(mpc)

            args = {};
            have_mpc = 0;
            if nargin > 0
                if isa(mpc, 'opf_model')
                    args = { mpc };
                elseif isstruct(mpc)
                    have_mpc = 1;
                end
            end

            %% call parent constructor
            om@opt_model(args{:});

            if have_mpc
                om.mpc = mpc;
            end

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if om.var is empty) in ADD_VAR(), DISPLAY() and
            %%              INIT_INDEXED_NAME(), after object construction,
            %%              but before object use.
        end

        function om = def_set_types(om)
            % Define set types ``var``, ``lin``, ``nle``, ``nli``, ``qdc``, ``nlc``, ``cost``.

            om.set_types = struct(...
                    'var',  'VARIABLES', ...
                    'lin',  'LINEAR CONSTRAINTS', ...
                    'nle',  'NONLIN EQ CONSTRAINTS', ...
                    'nli',  'NONLIN INEQ CONSTRAINTS', ...
                    'qdc',  'QUADRATIC COSTS', ...
                    'nlc',  'GEN NONLIN COSTS', ...
                    'cost', 'LEGACY COSTS'  ...
                );
        end

        function om = init_set_types(om)
            % Initialize data structures for each set type.

            %% call parent to create base data structures for each type
            init_set_types@opt_model(om);

            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            om.cost.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
            om.cost.params = [];
        end
    end     %% methods
end         %% classdef
