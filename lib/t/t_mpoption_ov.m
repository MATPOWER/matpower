function ov = t_mpoption_ov()
% t_mpoption_ov - Example of option overrides from file.

%   MATPOWER
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

ov = struct('verbose', 2, 'model', 'DC', 'opf', struct('dc', struct('solver', 'CPLEX')));
