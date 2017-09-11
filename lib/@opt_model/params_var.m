function [v0, vl, vu, vt] = params_var(om, name, idx)
%PARAMS_VAR  Returns initial value, lower bound and upper bound for opt variables.
%   [V0, VL, VU] = OM.PARAMS_VAR()
%   [V0, VL, VU] = OM.PARAMS_VAR(NAME)
%   [V0, VL, VU] = OM.PARAMS_VAR(NAME, IDX)
%   [V0, VL, VU, VT] = PARAMS_VAR(...)
%   Returns the initial value V0, lower bound VL and upper bound VU for
%   the full optimization variable vector, or for a specific named or named
%   and indexed variable set. Optionally also returns a corresponding char
%   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
%   integer and binary variables, respectively.
%
%   Examples:
%       [x, xmin, xmax] = om.params_var();
%       [Pg, Pmin, Pmax] = om.params_var('Pg');
%       [zij0, zijmin, zijmax, ztype] = om.params_var('z', {i, j});
%   
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
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
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields, updating only
    %% the subscripts before use
    sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
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
            % sc = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            sc(1).subs = name;
            sc(2).subs = idx;
            v0 = [ v0; subsref(om.var.data.v0, sc) ];
            vl = [ vl; subsref(om.var.data.vl, sc) ];
            vu = [ vu; subsref(om.var.data.vu, sc) ];
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % sn = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                sn = sc; sn(2).type = '()';
                N = subsref(om.var.idx.N, sn);
                vt0 = subsref(om.var.data.vt, sc);
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
            % sc = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
            v0 = subsref(om.var.data.v0, sc);
            vl = subsref(om.var.data.vl, sc);
            vu = subsref(om.var.data.vu, sc);
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % sn = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                sn = sc; sn(2).type = '()';
                N = subsref(om.var.idx.N, sn);
                vt0 = subsref(om.var.data.vt, sc);
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
