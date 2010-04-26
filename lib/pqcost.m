function [pcost, qcost] = pqcost(gencost, ng, on)
%PQCOST  Splits the gencost variable into two pieces if costs are given for Qg.
%   [PCOST, QCOST] = PQCOST(GENCOST, NG, ON) checks whether GENCOST has
%   cost information for reactive power generation (rows ng+1 to 2*ng).
%   If so, it returns the first NG rows in PCOST and the last NG rows in
%   QCOST. Otherwise, leaves QCOST empty. Also does some error checking.
%   If ON is specified (list of indices of generators which are on line)
%   it only returns the rows corresponding to these generators.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

if nargin < 3
    on = (1:ng)';
end

if size(gencost, 1) == ng
    pcost = gencost(on, :);
    qcost = [];
elseif size(gencost, 1) == 2 * ng
    pcost = gencost(on, :);
    qcost = gencost(on+ng, :);
else
    error('pqcost: gencost has wrong number of rows');
end
