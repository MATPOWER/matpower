function label = describe_idx(om, idx_type, idxs)
%DESCRIBE_IDX  Identifies variable, constraint and cost row indexes.
%   LABEL = DESCRIBE_IDX(OM, IDX_TYPE, IDXS)
%
%	Returns strings describing (name and index) the variable, constraint
%   or cost row that corresponds to the indices in IDXS. IDX_TYPE must be
%   one of the following: 'var', 'lin', 'nln', or 'cost', corresponding
%   to indices for variables, linear constraints, non-linear constraints
%	and cost rows, respectively. The return value is a string if IDXS is
%	a scalar, otherwise it is a cell array of strings of the same
% 	dimension as IDXS.
%
%   Examples:
%       label = describe_idx(om, 'var', 87));
%       labels = describe_idx(om, 'lin', [38; 49; 93]));
%   
%   See also OPT_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2012 by Power System Engineering Research Center (PSERC)
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

label = cell(size(idxs));		%% pre-allocate return cell array
for i = 1:length(idxs(:))
	ii = idxs(i);
	if ii > om.(idx_type).N
		error('@opt_model/describe_idx: index exceeds maximum %s index (%d)', idx_type, om.(idx_type).N);
	end
	if ii < 1
		error('@opt_model/describe_idx: index must be positive');
	end
	for k = om.(idx_type).NS:-1:1
		name = om.(idx_type).order(k).name;
		idx = om.(idx_type).order(k).idx;
		if isempty(idx)
			if ii >= om.(idx_type).idx.i1.(name)
				label{i} = sprintf('%s(%d)', name, ii - om.(idx_type).idx.i1.(name) + 1);
				break;
			end
		else
			s = substruct('.', name, '()', idx);
			if ii >= subsref(om.(idx_type).idx.i1, s)
				idxstr = sprintf('%d', idx{1});
				for j = 2:length(idx)
					idxstr = sprintf('%s,%d', idxstr, idx{j});
				end
				label{i} = sprintf('%s(%s)(%d)', name, idxstr, ...
							ii - subsref(om.(idx_type).idx.i1, s) + 1);
				break;
			end
		end
	end
end
if isscalar(idxs)				%% return scalar
	label = label{1};
end
