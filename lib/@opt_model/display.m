function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
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
