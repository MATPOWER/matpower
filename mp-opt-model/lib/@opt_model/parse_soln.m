function ps = parse_soln(om, stash)
% parse_soln - Parse solution vector and shadow prices by all named sets.
% ::
%
%   PS = OM.PARSE_SOLN()
%   OM.PARSE_SOLN(STASH)
%
%   For a solved model, OM.PARSE_SOLN() returns a struct of parsed
%   solution vector and shadow price values for each named set of
%   variables and constraints. The returned PS (parsed solution) struct
%   has the following format, where each of the terminal elements is a
%   struct with fields corresponding to the respective named sets:
%
%   Output:
%       PS
%           .var
%               .val
%               .mu_l
%               .mu_u
%           .lin
%               .mu_l
%               .mu_u
%           .nle
%               .lam
%           .nli
%               .mu
%
%   The value of each element in the returned struct can be obtained
%   via the GET_SOLN method as well, but using PARSE_SOLN is generally
%   more efficient if a complete set of values is needed.
%
%   If the optional STASH input argument is present and true, the fields
%   of the return struct are copied to OM.SOLN.
%
% See also get_soln.

%   MP-Opt-Model
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 2
    stash = false;
end

%% calls to substruct() are relatively expensive, so we pre-build the
%% structs for addressing cell and numeric array fields, updating only
%% the subscripts before use
sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
sn = sc; sn(2).type = '()';                         %% num array field

ps = struct('var', []);     %% parsed solution
s = om.soln;                %% aggregate solution

%% var
params = struct('src', s.x, 'dst', 'val');
if isfield(s.lambda, 'lower')
    params(end+1).src = s.lambda.lower;
    params(end  ).dst = 'mu_l';
end
if isfield(s.lambda, 'upper')
    params(end+1).src = s.lambda.upper;
    params(end  ).dst = 'mu_u';
end
ps.var = parse_soln_fields(om, 'var', params, sn, sc);

%% lin
if om.getN('lin')
    if isfield(s.lambda, 'mu_l')
        if isfield(s.lambda, 'mu_u')
            params = struct('src', {s.lambda.mu_l, s.lambda.mu_u}, ...
                            'dst', {'mu_l', 'mu_u'});
        else
            params = struct('src', s.lambda.mu_l, 'dst', 'mu_l');
        end
    else
        if isfield(s.lambda, 'mu_u')
            params = struct('src', s.lambda.mu_u, 'dst', 'mu_u');
        else
            params = [];
        end
    end
    if ~isempty(params)
        ps.lin = parse_soln_fields(om, 'lin', params, sn, sc);
    end
end

%% nle
if om.getN('nle') && isfield(s.lambda, 'eqnonlin')
    params = struct('src', s.lambda.eqnonlin, ...
                    'dst',     'lam'  );
    ps.nle = parse_soln_fields(om, 'nle', params, sn, sc);
end

%% nli
if om.getN('nli') && isfield(s.lambda, 'ineqnonlin')
    params = struct('src', s.lambda.ineqnonlin, ...
                    'dst',      'mu' );
    ps.nli = parse_soln_fields(om, 'nli', params, sn, sc);
end

%% stash the result, if requested
if stash
    om.soln = nested_struct_copy(om.soln, ps, struct('copy_mode', '='));
end



function psf = parse_soln_fields(om, ff, params, sn, sc)
om_ff = om.(ff);
psf = struct();     %% parsed solution field

np = length(params);
have_param = zeros(np, 1);
for j = 1:np
    have_param(j) = ~isempty(params(j).src);
end
for k = 1:om_ff.NS
    name = om_ff.order(k).name;
    idx = om_ff.order(k).idx;
    if isempty(idx)
        N = om_ff.idx.N.(name);
    else
        sn(1).subs = name;
        sn(2).subs = idx;
        N = subsref(om_ff.idx.N, sn);
        need_init = all([idx{:}] == 1);
    end
    if N
        for j = 1:np
            if have_param(j)    %% parameter is available
                dname = params(j).dst;  %% destination field name
                if isempty(idx)
                    i1 = om_ff.idx.i1.(name);
                    iN = om_ff.idx.iN.(name);
                    psf.(dname).(name)  = params(j).src(i1:iN);
                else
                    if need_init
                        switch ff
                            case 'var'
                                psf.(dname).(name) = cell(size(om_ff.data.v0.(name)));
                            case 'lin'
                                psf.(dname).(name) = cell(size(om_ff.data.A.(name)));
                            case 'qdc'
                                psf.(dname).(name) = cell(size(om_ff.data.c.(name)));
                            otherwise
                                psf.(dname).(name) = cell(size(om_ff.data.fcn.(name)));
                        end
                    end
                    i1 = subsref(om_ff.idx.i1, sn);    %% starting row index
                    iN = subsref(om_ff.idx.iN, sn);    %% ending row index
                    sc(1).subs = name;
                    sc(2).subs = idx;
                    psf.(dname) = ...
                        subsasgn(psf.(dname), sc, params(j).src(i1:iN));
                end
            end
        end     %% for j
    end
end     %% for k
