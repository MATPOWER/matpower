function ov = t_mpoption_ov()
%T_MPOPTION_OV  Example of option overrides from file.

%   MATPOWER
%   Copyright (c) 2013-2015 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

ov = struct('verbose', 2, 'model', 'DC', 'opf', struct('dc', struct('solver', 'CPLEX')));
