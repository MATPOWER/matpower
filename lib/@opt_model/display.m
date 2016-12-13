function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also OPT_MODEL.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
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
    fprintf('%10d = var.NS%29d = var.N\n', om.var.NS, om.var.N);
%     for k = 1:om.var.NS
%         name = om.var.order(k).name;
%         if isscalar(idx.N.(name))
%             fprintf('%10d:%19s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
%         else
%             d = size(idx.N.(name));     %% dimensions of this named set
%             nn = prod(d);               %% total number of blocks in named set
%             temp = cell(size(d));
%             [temp{end:-1:1}] = ind2sub(d(end:-1:1), (1:nn)');
%             ss = num2cell([temp{:}]);   %% table of indices
%             str = '%d'; for m = 2:length(d), str = [str ',%d']; end
%             for i = 1:nn
%                 s = substruct('.', name, '()', {ss{i,:}});
%                 nname = sprintf(['%s(' str, ')'], name, ss{i,:});
%                 fprintf('%10d:%19s %8d %8d %8d\n', k, nname, ...
%                     subsref(idx.i1, s), subsref(idx.iN, s), subsref(idx.N, s));
%             end
%         end
%     end
%     fprintf('%10s%38s\n', sprintf('var.NS = %d', om.var.NS), sprintf('var.N = %d', om.var.N));
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'VARIABLES');
end
if om.nln.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'NON-LINEAR CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '======================', '------', '-----', '-----', '------');
    idx = om.nln.idx;
    for k = 1:om.nln.NS
        name = om.nln.order(k).name;
        fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
    end
    fprintf('%10d = nln.NS%29d = nln.N\n', om.nln.NS, om.nln.N);
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'NON-LINEAR CONSTRAINTS');
end
if om.lin.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'LINEAR CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '==================', '------', '-----', '-----', '------');
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
    fprintf('%10d = lin.NS%29d = lin.N\n', om.lin.NS, om.lin.N);
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'LINEAR CONSTRAINTS');
end
if om.cost.NS
    fprintf('\n%-17s %12s %8s %8s %8s\n', 'COSTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-17s %12s %8s %8s %8s\n', '=====', '------', '-----', '-----', '------');
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
    fprintf('%10d = cost.NS%28d = cost.N\n', om.cost.NS, om.cost.N);
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'COSTS');
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
