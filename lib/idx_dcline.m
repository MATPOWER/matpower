function c = idx_dcline
%IDX_DCLINE   Defines constants for named column indices to dcline matrix.
%   Example:
%
%   c = idx_dcline;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    mpc.dcline(4, c.BR_STATUS) = 0;        % take branch 4 out of service
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-17 must be included in input matrix (in case file)
%    1  F_BUS     f, "from" bus number
%    2  T_BUS     t,  "to"  bus number
%    3  BR_STATUS initial branch status, 1 - in service, 0 - out of service
%    4  PF        MW flow at "from" bus ("from" -> "to")
%    5  PT        MW flow at  "to"  bus ("from" -> "to")
%    6  QF        MVAr injection at "from" bus ("from" -> "to")
%    7  QT        MVAr injection at  "to"  bus ("from" -> "to")
%    8  VF        voltage setpoint at "from" bus (p.u.)
%    9  VT        voltage setpoint at  "to"  bus (p.u.)
%   10  PMIN      lower limit on PF (MW flow at "from" end)
%   11  PMAX      upper limit on PF (MW flow at "from" end)
%   12  QMINF     lower limit on MVAr injection at "from" bus
%   13  QMAXF     upper limit on MVAr injection at "from" bus
%   14  QMINT     lower limit on MVAr injection at  "to"  bus
%   15  QMAXT     upper limit on MVAr injection at  "to"  bus
%   16  LOSS0     constant term of linear loss function (MW)
%   17  LOSS1     linear term of linear loss function (MW/MW)
%                 (loss = LOSS0 + LOSS1 * PF)
%
%   columns 18-23 are added to matrix after OPF solution
%   they are typically not present in the input matrix
%                 (assume OPF objective function has units, u)
%   18  MU_PMIN   Kuhn-Tucker multiplier on lower flow lim at "from" bus (u/MW)
%   19  MU_PMAX   Kuhn-Tucker multiplier on upper flow lim at "from" bus (u/MW)
%   20  MU_QMINF  Kuhn-Tucker multiplier on lower VAr lim at "from" bus (u/MVAr)
%   21  MU_QMAXF  Kuhn-Tucker multiplier on upper VAr lim at "from" bus (u/MVAr)
%   22  MU_QMINT  Kuhn-Tucker multiplier on lower VAr lim at  "to"  bus (u/MVAr)
%   23  MU_QMAXT  Kuhn-Tucker multiplier on upper VAr lim at  "to"  bus (u/MVAr)
%
%   See also TOGGLE_DCLINE.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2011 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define the indices
c = struct( ...
    'F_BUS',     1, ... %% f, "from" bus number
    'T_BUS',     2, ... %% t,  "to"  bus number
    'BR_STATUS', 3, ... %% initial branch status, 1 - in service, 0 - out of service
    'PF',        4, ... %% MW flow at "from" bus ("from" -> "to")
    'PT',        5, ... %% MW flow at  "to"  bus ("from" -> "to")
    'QF',        6, ... %% MVAr injection at "from" bus ("from" -> "to")
    'QT',        7, ... %% MVAr injection at  "to"  bus ("from" -> "to")
    'VF',        8, ... %% voltage setpoint at "from" bus (p.u.)
    'VT',        9, ... %% voltage setpoint at  "to"  bus (p.u.)
    'PMIN',     10, ... %% lower limit on PF (MW flow at "from" end)
    'PMAX',     11, ... %% upper limit on PF (MW flow at "from" end)
    'QMINF',    12, ... %% lower limit on MVAr injection at "from" bus
    'QMAXF',    13, ... %% upper limit on MVAr injection at "from" bus
    'QMINT',    14, ... %% lower limit on MVAr injection at  "to"  bus
    'QMAXT',    15, ... %% upper limit on MVAr injection at  "to"  bus
    'LOSS0',    16, ... %% constant term of linear loss function (MW)
    'LOSS1',    17, ... %% linear term of linear loss function (MW)
    'MU_PMIN',  18, ... %% Kuhn-Tucker multiplier on lower flow lim at "from" bus (u/MW)
    'MU_PMAX',  19, ... %% Kuhn-Tucker multiplier on upper flow lim at "from" bus (u/MW)
    'MU_QMINF', 20, ... %% Kuhn-Tucker multiplier on lower VAr lim at "from" bus (u/MVAr)
    'MU_QMAXF', 21, ... %% Kuhn-Tucker multiplier on upper VAr lim at "from" bus (u/MVAr)
    'MU_QMINT', 22, ... %% Kuhn-Tucker multiplier on lower VAr lim at  "to"  bus (u/MVAr)
    'MU_QMAXT', 23  );  %% Kuhn-Tucker multiplier on upper VAr lim at  "to"  bus (u/MVAr)
