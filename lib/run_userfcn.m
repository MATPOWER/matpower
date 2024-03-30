function rv = run_userfcn(userfcn, stage, varargin)
% run_userfcn - Runs the userfcn callbacks for a given stage.
% ::
%
%   RV = RUN_USERFCN(USERFCN, STAGE, VARARGIN)
%
%   USERFCN : the 'userfcn' field of mpc, populated by ADD_USERFCN
%   STAGE   : the name of the callback stage being executed
%   (additional arguments) some stages require additional arguments.
%
%   Example:
%       mpc = om.get_mpc();
%       om = run_userfcn(mpc.userfcn, 'formulation', om);
%
% See also add_userfcn, remove_userfcn, toggle_reserves, toggle_iflims,
% runopf_w_res.

%   MATPOWER
%   Copyright (c) 2009-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

rv = varargin{1};
if ~isempty(userfcn) && isfield(userfcn, stage)
    for k = 1:length(userfcn.(stage))
        if isfield(userfcn.(stage)(k), 'args')
            args = userfcn.(stage)(k).args;
        else
            args = [];
        end
        rv = feval(userfcn.(stage)(k).fcn, rv, varargin{2:end}, args);
        % mpc     = userfcn_*_ext2int(mpc, mpopt, args);
        % om      = userfcn_*_formulation(om, mpopt, args);
        % results = userfcn_*_int2ext(results, mpopt, args);
        % results = userfcn_*_printpf(results, fd, mpopt, args);
        % mpc     = userfcn_*_savecase(mpc, fd, prefix, args);
    end
end
