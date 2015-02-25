function [v0, vl, vu] = getv(om, name, idx)
%GETV  Returns initial value, lower bound and upper bound for opt variables.
%   [V0, VL, VU] = GETV(OM)
%   [V0, VL, VU] = GETV(OM, NAME)
%   [V0, VL, VU] = GETV(OM, NAME, IDX)
%   Returns the lower bound and upper bound for the full optimization
%   variable vector, or for a specific named or name and indexed
%   variable set.
%
%   Examples:
%       [x, xmin, xmax] = getv(om);
%       [Pg, Pmin, Pmax] = getv(om, 'Pg');
%       [zij0, zijmin, zijmax] = getv(om, 'z', {i, j});
%   
%   See also OPT_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://matpower.org/ for more info.

if nargin < 2
    v0 = []; vl = []; vu = [];
    for k = 1:om.var.NS
        name = om.var.order(k).name;
        idx = om.var.order(k).idx;
        if isempty(idx)
            v0 = [ v0; om.var.data.v0.(name) ];
            vl = [ vl; om.var.data.vl.(name) ];
            vu = [ vu; om.var.data.vu.(name) ];
        else
            s = substruct('.', name, '{}', idx);
            v0 = [ v0; subsref(om.var.data.v0, s) ];
            vl = [ vl; subsref(om.var.data.vl, s) ];
            vu = [ vu; subsref(om.var.data.vu, s) ];
        end
    end
else
    if isfield(om.var.idx.N, name)
        if nargin < 3
            v0 = om.var.data.v0.(name);
            vl = om.var.data.vl.(name);
            vu = om.var.data.vu.(name);
        else
            s1 = substruct('.', name, '{}', idx);
            v0 = subsref(om.var.data.v0, s1);
            vl = subsref(om.var.data.vl, s1);
            vu = subsref(om.var.data.vu, s1);
        end
    else
        v0 = [];
        vl = [];
        vu = [];
    end
end
