function mpc = add_userfcn(mpc, stage, fcn, args, allow_multiple)
%ADD_USERFCN   Appends a userfcn to the list to be called for a case.
%
%   mpc = add_userfcn(mpc, stage, fcn)
%   mpc = add_userfcn(mpc, stage, fcn, args)
%   mpc = add_userfcn(mpc, stage, fcn, args, allow_multiple)
%
%   A userfcn is a callback function that can be called automatically by
%   MATPOWER at one of various stages in a simulation.
%
%   mpc   : the case struct
%   stage : the name of the stage at which this function should be
%           called: ext2int, formulation, int2ext, printpf
%   fcn   : the name of the userfcn
%   args  : (optional) the value to be passed as an argument to the
%           userfcn (typically a struct)
%   allow_multiple : (optional) if TRUE, allows the same function to
%          be added more than once.
%
%   Currently there are 5 different callback stages defined. Each stage has
%   a name, and by convention, the name of a user-defined callback function
%   ends with the name of the stage. The following is a description of each
%   stage, when it is called and the input and output arguments (which vary
%   depending on the stage). The reserves example (see 'help runopf_w_res')
%   is used to illustrate how these callback userfcn's might be used.
%
%   1. ext2int
%
%   Called from ext2int() immediately after the case is converted from
%   external to internal indexing. Inputs are a MATPOWER case struct (mpc),
%   freshly converted to internal indexing and any (optional) args value
%   supplied to add_userfcn. Output is the (presumably updated) mpc. This is
%   typically used to reorder any input arguments that may be needed in
%   internal ordering by the formulation stage.
%
%   E.g. mpc = userfcn_reserves_ext2int(mpc, args)
%
%   2. formulation
%
%   Called from opf() after the OPF Model (OM) object has been initialized
%   with the standard OPF formulation, but before calling the solver. Inputs
%   are the OM object and any (optional) args supplied to add_userfcn.
%   Output is the om object. This is the ideal place to add any additional
%   vars, constraints or costs to the OPF formulation.
%
%   E.g. om = userfcn_reserves_formulation(om, args)
%
%   3. int2ext
%
%   Called from int2ext() immediately before the resulting case is converted
%   from internal back to external indexing. Inputs are the results struct
%   and any (optional) args supplied via add_userfcn. Output is the results
%   struct. This is typically used to convert any results to external
%   indexing and populate any corresponding fields in the results struct.
%
%   E.g. results = userfcn_reserves_int2ext(results, args)
%
%   4. printpf
%
%   Called from printpf() after the pretty-printing of the standard OPF
%   output. Inputs are the results struct, the file descriptor to write to,
%   a MATPOWER options vector, and any (optional) args supplied via
%   add_userfcn. Output is the results struct. This is typically used for
%   any additional pretty-printing of results.
%
%   E.g. results = userfcn_reserves_printpf(results, fd, mpopt, args)
%
%   5. savecase
%
%   Called from savecase() when saving a case struct to an m-file after
%   printing all of the other data to the file. Inputs are the case struct,
%   the file descriptor to write to, the variable prefix (typically 'mpc.')
%   and any (optional) args supplied via add_userfcn. Output is the case
%   struct. This is typically used to write any non-standard case struct
%   fields to the case file.
%
%   E.g. mpc = userfcn_reserves_printpf(mpc, fd, prefix, args)

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
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
            for k = 1:n-1
                if strcmp(mpc.userfcn.(stage)(k).fcn, fcn)
                    error('add_userfcn: the function ''%s'' has already been added', fcn);
                end
            end
        end
    end
end

mpc.userfcn.(stage)(n).fcn = fcn;
if ~isempty(args)
    mpc.userfcn.(stage)(n).args = args;
end

return;
