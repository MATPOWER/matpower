function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

display@opt_model(om)

if om.cost.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'LEGACY COSTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '============', '------', '-----', '-----', '------');
    idx = om.cost.idx;
    for k = 1:om.cost.NS
        name = om.cost.order(k).name;
        if isempty(om.cost.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.cost.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = cost.NS%28d = cost.N\n\n', om.cost.NS, om.cost.N);
else
    fprintf('%s  :  <none>\n', 'LEGACY COSTS');
end

fprintf('  mpc = ');
if ~isempty(fieldnames(om.mpc))
    fprintf('\n');
end
if have_fcn('octave')
    fprintf('    <scalar struct>\n');
else
    display(om.mpc);
end
