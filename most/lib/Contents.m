% MOST
%   Version 1.0.1       30-Oct-2018
%
%   The MATPOWER Optimal Scheduling Tool (MOST) is framework for solving
%   generalized steady-state electric power scheduling problems.
%
%   MOST can be used to solve problems as simple as a deterministic,
%   single period economic dispatch problem with no transmission
%   constraints or as complex as a stochastic, security-constrained,
%   combined unit-commitment and multiperiod optimal power flow problem
%   with locational contingency and load-following reserves, ramping
%   costs and constraints, deferrable demands, lossy storage resources
%   and uncertain renewable generation.
%
%   While the problem formulation is general and incorporates a full
%   nonlinear AC network model, the current implementation is limited
%   to DC power flow modeling of the network. Some work has been done
%   on an AC implementation, but it is not yet ready for release.
%
%   The primary developers of MOST are Carlos E. Murillo-Sanchez and
%   Ray D. Zimmerman, with significant contributions from
%   Daniel Munoz-Alvarez and Alberto J. Lamadrid. It is built on top of
%   MATPOWER, a package of MATLAB/Octave M-files for solving power flow
%   and optimal power flow problems.
%
%   The latest version of MOST can be found on GitHub at:
%
%       https://github.com/MATPOWER/most
%
%   MOST is covered by the 3-clause BSD License (see LICENSE for details).

%   MOST
%   Copyright (c) 2010-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.
