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
