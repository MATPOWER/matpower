function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also OPF_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2010 by Power System Engineering Research Center (PSERC)
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
    fprintf('\n%-22s %5s %8s %8s %8s\n', 'VARIABLES', 'name', 'i1', 'iN', 'N');
    fprintf('%-22s %5s %8s %8s %8s\n', '=========', '------', '-----', '-----', '------');
    for k = 1:om.var.NS
        name = om.var.order{k};
        idx = om.var.idx;
        fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
    end
    fprintf('%15s%31s\n', sprintf('var.NS = %d', om.var.NS), sprintf('var.N = %d', om.var.N));
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'VARIABLES');
end
if om.nln.NS
    fprintf('\n%-22s %5s %8s %8s %8s\n', 'NON-LINEAR CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-22s %5s %8s %8s %8s\n', '======================', '------', '-----', '-----', '------');
    for k = 1:om.nln.NS
        name = om.nln.order{k};
        idx = om.nln.idx;
        fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
    end
    fprintf('%15s%31s\n', sprintf('nln.NS = %d', om.nln.NS), sprintf('nln.N = %d', om.nln.N));
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'NON-LINEAR CONSTRAINTS');
end
if om.lin.NS
    fprintf('\n%-22s %5s %8s %8s %8s\n', 'LINEAR CONSTRAINTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-22s %5s %8s %8s %8s\n', '==================', '------', '-----', '-----', '------');
    for k = 1:om.lin.NS
        name = om.lin.order{k};
        idx = om.lin.idx;
        fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
    end
    fprintf('%15s%31s\n', sprintf('lin.NS = %d', om.lin.NS), sprintf('lin.N = %d', om.lin.N));
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'LINEAR CONSTRAINTS');
end
if om.cost.NS
    fprintf('\n%-22s %5s %8s %8s %8s\n', 'COSTS', 'name', 'i1', 'iN', 'N');
    fprintf('%-22s %5s %8s %8s %8s\n', '=====', '------', '-----', '-----', '------');
    for k = 1:om.cost.NS
        name = om.cost.order{k};
        idx = om.cost.idx;
        fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
    end
    fprintf('%15s%31s\n', sprintf('cost.NS = %d', om.cost.NS), sprintf('cost.N = %d', om.cost.N));
    fprintf('\n');
else
    fprintf('%s  :  <none>\n', 'COSTS');
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

fprintf('  userdata = ');
if ~isempty(fieldnames(om.userdata))
    fprintf('\n');
end
if have_fcn('octave')
    fprintf('    <scalar struct>\n');
else
    display(om.userdata);
end
