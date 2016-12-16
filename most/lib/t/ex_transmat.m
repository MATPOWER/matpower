function transmat = ex_transmat(nt)
%EX_TRANSMAT  Example transition probability matrix definition.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

transmat = cell(1, nt);
T = [ 0.158655253931457; 0.682689492137086; 0.158655253931457 ];
[transmat{:}] = deal(T * ones(1,3));
transmat{1} = T;
