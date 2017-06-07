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
%   This class is a sub-class of OPT_MODEL and simply adds the 'mpc'
%   field for storing the MATPOWER case struct used to build the object
%   along with the get_mpc() method.
%
%   The following is the structure of the data in the OPF model object.
%
%   om
%       <opt_model fields> - see OPT_MODEL for details
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
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

    properties
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
