function str = display(om, name, N)
%DISPLAY  Displays the object
%

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.


fprintf('%16s: %g\n', 'var.N', om.var.N);
fprintf('%16s: %g\n', 'nln.N', om.nln.N);
fprintf('%16s: %g\n', 'lin.N', om.lin.N);
fprintf('%16s: %g\n', 'var.NS', om.var.NS);
fprintf('%16s: %g\n', 'nln.NS', om.nln.NS);
fprintf('%16s: %g\n', 'lin.NS', om.lin.NS);
if om.var.NS
	fprintf('\n%16s: %10s %8s %8s %8s\n', 'vars', 'name', 'i1', 'iN', 'N');
	fprintf('%29s %8s %8s %8s\n', '------', '-----', '-----', '------');
	for k = 1:om.var.NS
		name = om.var.order{k};
		idx = om.var.idx;
		fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
	end
else
	fprintf('%16s: []\n', 'varsets');
end
if om.nln.NS
	fprintf('\n%16s: %10s %8s %8s %8s\n', 'nln cons', 'name', 'i1', 'iN', 'N');
	fprintf('%29s %8s %8s %8s\n', '------', '-----', '-----', '------');
	for k = 1:om.nln.NS
		name = om.nln.order{k};
		idx = om.nln.idx;
		fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
	end
else
	fprintf('%16s: []\n', 'nlconsets');
end
if om.lin.NS
	fprintf('\n%16s: %10s %8s %8s %8s\n', 'lin cons', 'name', 'i1', 'iN', 'N');
	fprintf('%29s %8s %8s %8s\n', '------', '-----', '-----', '------');
	for k = 1:om.lin.NS
		name = om.lin.order{k};
		idx = om.lin.idx;
		fprintf('%15d:%12s %8d %8d %8d\n', k, name, idx.i1.(name), idx.iN.(name), idx.N.(name));
	end
else
	fprintf('%16s: []\n', 'lnconsets');
end

if ~isempty(om.mpc)
	mpc = om.mpc
else
	fprintf('%16s: []\n', 'mpc');
end

return;
