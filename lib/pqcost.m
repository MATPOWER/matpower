function [pcost, qcost] = pqcost(gencost, ng, on)
%PQCOST  Splits the gencost variable into two pieces if costs are given for Qg.
%   [pcost, qcost] = pqcost(gencost, ng, on) checks whether gencost has cost
%   information for reactive power generation (rows ng+1 to 2*ng). If so it
%   returns the first ng rows in pcost and the last ng rows in qcost.
%   Otherwise, leaves qcost empty. Also does some error checking.
%   If on is specified (list of indices of generators which are on line)
%   it only returns the rows corresponding to these generators.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 3
    on = [1:ng]';
end

if size(gencost, 1) == ng
    pcost = gencost(on, :);
    qcost = [];
elseif size(gencost, 1) == 2 * ng
    pcost = gencost(on, :);
    qcost = gencost(on+ng, :);
else
    error('gencost has wrong number of rows');
end

return;
