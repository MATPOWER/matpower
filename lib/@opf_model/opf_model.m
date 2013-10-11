function om = opf_model(mpc)
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
%   Each field of .idx or .data is a struct whose field names are the names
%   of the corresponding blocks of vars, constraints or costs (found in
%   order in the corresponding .order field). The description next to these
%   fields gives the meaning of the value for each named sub-field.
%   E.g. om.var.data.v0.Pg contains a vector of initial values for the 'Pg'
%   block of variables.
%
%   om
%       .opt_model  - the corresponding OPT_MODEL object
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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin == 0
	es = struct();
    s = struct('mpc', es);
    om = opt_model;
    om = class(s, 'opf_model', om);
else
    if isa(mpc,'opf_model') 
        om = mpc;
    else
        s = struct('mpc', mpc);
        om = opt_model;
        om = class(s, 'opf_model', om);
    end
end
