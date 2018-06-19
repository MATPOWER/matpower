function newval = e2i_data(mpc, val, ordering, dim)
%E2I_DATA   Converts data from external to internal indexing.
%
%   VAL = E2I_DATA(MPC, VAL, ORDERING)
%   VAL = E2I_DATA(MPC, VAL, ORDERING, DIM)
%
%   When given a case struct that has already been converted to
%   internal indexing, this function can be used to convert other data
%   structures as well by passing in 2 or 3 extra parameters in
%   addition to the case struct. If the value passed in the 2nd
%   argument is a column vector or cell array, it will be converted
%   according to the ORDERING specified by the 3rd argument (described
%   below). If VAL is an n-dimensional matrix or cell array, then the
%   optional 4th argument (DIM, default = 1) can be used to specify
%   which dimension to reorder. The return value in this case is the
%   value passed in, converted to internal indexing.
%
%   The 3rd argument, ORDERING, is used to indicate whether the data
%   corresponds to bus-, gen- or branch-ordered data. It can be one
%   of the following three strings: 'bus', 'gen' or 'branch'. For
%   data structures with multiple blocks of data, ordered by bus,
%   gen or branch, they can be converted with a single call by
%   specifying ORDERING as a cell array of strings.
%
%   Any extra elements, rows, columns, etc. beyond those indicated
%   in ORDERING, are not disturbed.
%
%   Examples:
%       A_int = e2i_data(mpc, A_ext, {'bus','bus','gen','gen'}, 2);
%
%       Converts an A matrix for user-supplied OPF constraints from
%       external to internal ordering, where the columns of the A
%       matrix correspond to bus voltage angles, then voltage
%       magnitudes, then generator real power injections and finally
%       generator reactive power injections.
%
%       gencost_int = e2i_data(mpc, gencost_ext, {'gen','gen'}, 1);
%
%       Converts a GENCOST matrix that has both real and reactive power
%       costs (in rows 1--ng and ng+1--2*ng, respectively).
%
%   See also I2E_DATA, E2I_FIELD, EXT2INT.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if ~isfield(mpc, 'order')
    error('e2i_data: mpc does not have the ''order'' field required to convert from external to internal numbering.');
end
o = mpc.order;
if o.state ~= 'i'
    error('e2i_data: mpc does not have internal ordering data available, call ext2int first');
end
if nargin < 4
    dim = 1;
end
if ischar(ordering)         %% single set
    if strcmp(ordering, 'gen')
        idx = o.(ordering).status.on(o.(ordering).i2e);
    else
        idx = o.(ordering).status.on;
    end
    newval = get_reorder(val, idx, dim);
else                        %% multiple sets
    b = 0;  %% base
    for k = 1:length(ordering)
        n = size(o.ext.(ordering{k}), 1);
        v = get_reorder(val, b+(1:n), dim);
        new_v{k} = e2i_data(mpc, v, ordering{k}, dim);
        b = b + n;
    end
    n = size(val, dim);
    if n > b                %% the rest
        v = get_reorder(val, b+1:n, dim);
        new_v{length(new_v)+1} = v;
    end
    newval = cat(dim, new_v{:});
end
