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
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
