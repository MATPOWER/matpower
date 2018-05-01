function mpc = add_userfcn(mpc, stage, fcn, args, allow_multiple)
%ADD_USERFCN   Appends a userfcn to the list to be called for a case.
%
%   MPC = ADD_USERFCN(MPC, STAGE, FCN)
%   MPC = ADD_USERFCN(MPC, STAGE, FCN, ARGS)
%   MPC = ADD_USERFCN(MPC, STAGE, FCN, ARGS, ALLOW_MULTIPLE)
%
%   A userfcn is a callback function that can be called automatically by
%   MATPOWER at one of various stages in a simulation.
%
%   MPC   : the case struct
%   STAGE : the name of the stage at which this function should be
%           called: ext2int, formulation, int2ext, printpf, savecase
%   FCN   : the name of the userfcn
%   ARGS  : (optional) the value to be passed as an argument to the
%           userfcn (typically a struct)
%   ALLOW_MULTIPLE : (optional) if TRUE, allows the same function to
%          be added more than once.
%
%   Currently there are 5 different callback stages defined. Each stage has
%   a name, and by convention, the name of a user-defined callback function
%   ends with the name of the stage. The following is a description of each
%   stage, when it is called and the input and output arguments which vary
%   depending on the stage. The reserves example (see RUNOPF_W_RES) is used
%   to illustrate how these callback userfcns might be used.
%
%   1. ext2int
%
%   Called from EXT2INT immediately after the case is converted from
%   external to internal indexing. Inputs are a MATPOWER case struct (MPC),
%   freshly converted to internal indexing and any (optional) ARGS value
%   supplied via ADD_USERFCN. Output is the (presumably updated) MPC. This is
%   typically used to reorder any input arguments that may be needed in
%   internal ordering by the formulation stage.
%
%   E.g. mpc = userfcn_reserves_ext2int(mpc, mpopt, args)
%
%   2. formulation
%
%   Called from OPF after the OPF Model (OM) object has been initialized
%   with the standard OPF formulation, but before calling the solver. Inputs
%   are the OM object and any (optional) ARGS supplied via ADD_USERFCN.
%   Output is the OM object. This is the ideal place to add any additional
%   vars, constraints or costs to the OPF formulation.
%
%   E.g. om = userfcn_reserves_formulation(om, mpopt, args)
%
%   3. int2ext
%
%   Called from INT2EXT immediately before the resulting case is converted
%   from internal back to external indexing. Inputs are the RESULTS struct
%   and any (optional) ARGS supplied via ADD_USERFCN. Output is the RESULTS
%   struct. This is typically used to convert any results to external
%   indexing and populate any corresponding fields in the RESULTS struct.
%
%   E.g. results = userfcn_reserves_int2ext(results, mpopt, args)
%
%   4. printpf
%
%   Called from PRINTPF after the pretty-printing of the standard OPF
%   output. Inputs are the RESULTS struct, the file descriptor to write to,
%   a MATPOWER options struct, and any (optional) ARGS supplied via
%   ADD_USERFCN. Output is the RESULTS struct. This is typically used for
%   any additional pretty-printing of results.
%
%   E.g. results = userfcn_reserves_printpf(results, fd, mpopt, args)
%
%   5. savecase
%
%   Called from SAVECASE when saving a case struct to an M-file after
%   printing all of the other data to the file. Inputs are the case struct,
%   the file descriptor to write to, the variable prefix (typically 'mpc.')
%   and any (optional) ARGS supplied via ADD_USERFCN. Output is the case
%   struct. This is typically used to write any non-standard case struct
%   fields to the case file.
%
%   E.g. mpc = userfcn_reserves_printpf(mpc, fd, prefix, args)
%
%   See also RUN_USERFCN, REMOVE_USERFCN, TOGGLE_RESERVES, TOGGLE_IFLIMS,
%   TOGGLE_DCLINE, TOGGLE_SOFTLIMS, and RUNOPF_W_RES.

%   MATPOWER
%   Copyright (c) 2009-2018, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 5
    allow_multiple = 0;
    if nargin < 4
        args = [];
    end
end
switch stage
    case {'ext2int', 'formulation', 'int2ext', 'printpf', 'savecase'}
        %% ok
    otherwise
        error('add_userfcn : ''%s'' is not the name of a valid callback stage\n', stage);
end

n = 1;
if isfield(mpc, 'userfcn')
    if isfield(mpc.userfcn, stage)
        n = length(mpc.userfcn.(stage)) + 1;
        if ~allow_multiple
            if have_fcn('octave')
                fcn_info = functions(fcn);
                for k = 1:n-1
                    cb_info = functions(mpc.userfcn.(stage)(k).fcn);
                    if strcmp(cb_info.function, fcn_info.function)
                        error('add_userfcn: the function ''%s'' has already been added', func2str(fcn));
                    end
                end
            else
                for k = 1:n-1
                    if isequal(mpc.userfcn.(stage)(k).fcn, fcn)
                        error('add_userfcn: the function ''%s'' has already been added', func2str(fcn));
                    end
                end
            end
        end
    end
end

mpc.userfcn.(stage)(n).fcn = fcn;
if ~isempty(args)
    mpc.userfcn.(stage)(n).args = args;
end
