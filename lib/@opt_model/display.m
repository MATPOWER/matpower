function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if om.var.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'VARIABLES', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '=========', '------', '-----', '-----', '------');
    idx = om.var.idx;
    for k = 1:om.var.NS
        name = om.var.order(k).name;
        if isempty(om.var.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.var.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = var.NS%29d = var.N\n\n', om.var.NS, om.var.N);
else
    fprintf('%s  :  <none>\n', 'VARIABLES');
end
if om.nle.NS
    fprintf('\n%-21s %8s %8s %8s %8s\n', 'NONLIN EQ CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-21s %8s %8s %8s %8s\n', '=====================', '------', '-----', '-----', '------');
    idx = om.nle.idx;
    for k = 1:om.nle.NS
        name = om.nle.order(k).name;
        if isempty(om.nle.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.nle.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = nle.NS%29d = nle.N\n\n', om.nle.NS, om.nle.N);
else
    fprintf('%s  :  <none>\n', 'NONLIN EQ CONSTRAINTS');
end
if om.nli.NS
    fprintf('\n%-23s %6s %8s %8s %8s\n', 'NONLIN INEQ CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-23s %6s %8s %8s %8s\n', '=======================', '------', '-----', '-----', '------');
    idx = om.nli.idx;
    for k = 1:om.nli.NS
        name = om.nli.order(k).name;
        if isempty(om.nli.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.nli.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = nli.NS%29d = nli.N\n\n', om.nli.NS, om.nli.N);
else
    fprintf('%s  :  <none>\n', 'NONLIN INEQ CONSTRAINTS');
end
if om.lin.NS
    fprintf('\n%-18s %11s %8s %8s %8s\n', 'LINEAR CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-18s %11s %8s %8s %8s\n', '==================', '------', '-----', '-----', '------');
    idx = om.lin.idx;
    for k = 1:om.lin.NS
        name = om.lin.order(k).name;
        if isempty(om.lin.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.lin.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = lin.NS%29d = lin.N\n\n', om.lin.NS, om.lin.N);
else
    fprintf('%s  :  <none>\n', 'LINEAR CONSTRAINTS');
end
if om.qdc.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'QUADRATIC COSTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '===============', '------', '-----', '-----', '------');
    idx = om.qdc.idx;
    for k = 1:om.qdc.NS
        name = om.qdc.order(k).name;
        if isempty(om.qdc.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.qdc.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = qdc.NS%28d = qdc.N\n\n', om.qdc.NS, om.qdc.N);
else
    fprintf('%s  :  <none>\n', 'QUADRATIC COSTS');
end
if om.nlc.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'GEN NONLIN COSTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '================', '------', '-----', '-----', '------');
    idx = om.nlc.idx;
    for k = 1:om.nlc.NS
        name = om.nlc.order(k).name;
        if isempty(om.nlc.order(k).idx)
            fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
        else
            vsidx = om.nlc.order(k).idx;
            str = '%d'; for m = 2:length(vsidx), str = [str ',%d']; end
            s = substruct('.', name, '()', vsidx);
            nname = sprintf(['%s(' str, ')'], name, vsidx{:});
            fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
                    subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
        end
    end
    fprintf('%10d = nlc.NS%28d = nlc.N\n\n', om.nlc.NS, om.nlc.N);
else
    fprintf('%s  :  <none>\n', 'GEN NONLIN COSTS');
end
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

fprintf('  userdata = ');
if ~isempty(fieldnames(om.userdata))
    fprintf('\n');
end
if have_fcn('octave')
    fprintf('    <scalar struct>\n');
else
    display(om.userdata);
end
