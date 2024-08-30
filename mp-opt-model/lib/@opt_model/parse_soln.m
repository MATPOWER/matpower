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

%% var
ps = struct('var', om.var.parse_soln(om.soln, stash));

%% lin
ps_lin = om.lin.parse_soln(om.soln, stash);
if ~isempty(ps_lin)
    ps.lin = ps_lin;
end

%% nle
ps_nle = om.nle.parse_soln(om.soln, true, stash);
if ~isempty(ps_nle)
    ps.nle = ps_nle;
end

%% nli
ps_nli = om.nli.parse_soln(om.soln, false, stash);
if ~isempty(ps_nli)
    ps.nli = ps_nli;
end

%%-----  DEPRECATED  -----
%% stash the result, if requested
if stash
    om.soln = nested_struct_copy(om.soln, ps, struct('copy_mode', '='));
end
