function [v0, vl, vu, vt] = getv(om, name, idx)
%GETV  Returns initial value, lower bound and upper bound for opt variables.
%   [V0, VL, VU] = GETV(OM)
%   [V0, VL, VU] = GETV(OM, NAME)
%   [V0, VL, VU] = GETV(OM, NAME, IDX)
%   [V0, VL, VU, VT] = GETV(...)
%   Returns the initial value V0, lower bound VL and upper bound VU for
%   the full optimization variable vector, or for a specific named or named
%   and indexed variable set. Optionally also returns a corresponding char
%   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
%   integer and binary variables, respectively.
%
%   Examples:
%       [x, xmin, xmax] = getv(om);
%       [Pg, Pmin, Pmax] = getv(om, 'Pg');
%       [zij0, zijmin, zijmax, ztype] = getv(om, 'z', {i, j});
%   
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargout > 3
    have_vt = 1;
else
    have_vt = 0;
end
if nargin < 2
    v0 = []; vl = []; vu = []; vt = char([]);
    s1 = struct('type', {'.', '{}'}, 'subs', {'', 1});
    for k = 1:om.var.NS
        name = om.var.order(k).name;
        idx = om.var.order(k).idx;
        if isempty(idx)
            v0 = [ v0; om.var.data.v0.(name) ];
            vl = [ vl; om.var.data.vl.(name) ];
            vu = [ vu; om.var.data.vu.(name) ];
            if have_vt
                N = om.var.idx.N.(name);
                vt0 = om.var.data.vt.(name);
                if isscalar(vt0) && N > 1 
                    vt = [ vt char(vt0 * ones(1, N)) ];
                else
                    vt = [ vt vt0 ];
                end
            end
        else
            % (calls to substruct() are relatively expensive ...
            % s1 = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            s1(1).subs = name;
            s1(2).subs = idx;
            v0 = [ v0; subsref(om.var.data.v0, s1) ];
            vl = [ vl; subsref(om.var.data.vl, s1) ];
            vu = [ vu; subsref(om.var.data.vu, s1) ];
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % s2 = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                s2 = s1;
                s2(2).type = '()';
                N = subsref(om.var.idx.N, s2);
                vt0 = subsref(om.var.data.vt, s1);
                if isscalar(vt0) && N > 1 
                    vt = [ vt char(vt0 * ones(1, N)) ];
                else
                    if ~isempty(vt0)
                        vt = [ vt vt0 ];
                    end
                end
            end
        end
    end
else
    if isfield(om.var.idx.N, name)
        if nargin < 3
            v0 = om.var.data.v0.(name);
            vl = om.var.data.vl.(name);
            vu = om.var.data.vu.(name);
            if have_vt
                N = om.var.idx.N.(name);
                vt0 = om.var.data.vt.(name);
                if isscalar(vt0) && N > 1 
                    vt = char(vt0 * ones(1, N));
                else
                    vt = vt0;
                end
            end
        else
            % (calls to substruct() are relatively expensive ...
            % s1 = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            s1 = struct('type', {'.', '{}'}, 'subs', {name, idx});
            v0 = subsref(om.var.data.v0, s1);
            vl = subsref(om.var.data.vl, s1);
            vu = subsref(om.var.data.vu, s1);
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % s2 = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                s2 = s1;
                s2(2).type = '()';
                N = subsref(om.var.idx.N, s2);
                vt0 = subsref(om.var.data.vt, s1);
                if isscalar(vt0) && N > 1 
                    vt = char(vt0 * ones(1, N));
                else
                    vt = vt0;
                end
            end
        end
    else
        v0 = [];
        vl = [];
        vu = [];
        if have_vt
            vt = [];
        end
    end
end
