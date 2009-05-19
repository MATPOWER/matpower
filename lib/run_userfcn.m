function rv = run_userfcn(userfcn, stage, varargin)
%RUN_USERFCN   Runs the userfcn callbacks for a given stage.
%
%   rv = run_userfcn(userfcn, stage, varargin)
%
%   userfcn : the 'userfcn' field of mpc, populated by add_userfcn()
%   stage   : the name of the callback stage begin executed
%   (additional arguments) some stages require additional arguments.
%
%   See 'help add_userfcn' for more details.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
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

return;
