function [pcost, qcost] = pqcost(gencost, ng, on)
%PQCOST  Splits the gencost variable into two pieces if costs are given for Qg.
%   [PCOST, QCOST] = PQCOST(GENCOST, NG, ON) checks whether GENCOST has
%   cost information for reactive power generation (rows ng+1 to 2*ng).
%   If so, it returns the first NG rows in PCOST and the last NG rows in
%   QCOST. Otherwise, leaves QCOST empty. Also does some error checking.
%   If ON is specified (list of indices of generators which are on line)
%   it only returns the rows corresponding to these generators.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

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
