function rv = run_userfcn(userfcn, stage, varargin)
%RUN_USERFCN   Runs the userfcn callbacks for a given stage.
%   RV = RUN_USERFCN(USERFCN, STAGE, VARARGIN)
%
%   USERFCN : the 'userfcn' field of mpc, populated by ADD_USERFCN
%   STAGE   : the name of the callback stage begin executed
%   (additional arguments) some stages require additional arguments.
%
%   Example:
%       mpc = get_mpc(om);
%       om = run_userfcn(mpc.userfcn, 'formulation', om);
%
%   See also ADD_USERFCN, REMOVE_USERFCN, TOGGLE_RESERVES, TOGGLE_IFLIMS,
%   RUNOPF_W_RES.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
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
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

rv = varargin{1};
if ~isempty(userfcn) && isfield(userfcn, stage)
    for k = 1:length(userfcn.(stage))
        if isfield(userfcn.(stage)(k), 'args')
            args = userfcn.(stage)(k).args;
        else
            args = [];
        end
        switch stage
            case {'ext2int', 'formulation', 'int2ext'}
                % mpc     = userfcn_*_ext2int(mpc, args);
                % om      = userfcn_*_formulation(om, args);
                % results = userfcn_*_int2ext(results, args);
                rv = feval(userfcn.(stage)(k).fcn, rv, args);
            case {'printpf', 'savecase'}
                % results = userfcn_*_printpf(results, fd, mpopt, args);
                % mpc     = userfcn_*_savecase(mpc, fd, prefix, args);
                rv = feval(userfcn.(stage)(k).fcn, rv, varargin{2}, varargin{3}, args);
        end
    end
end
