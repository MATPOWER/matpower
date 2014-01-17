function opt = mpoption_info_gurobi(selector)
%MPOPTION_INFO_GUROBI  Returns MATPOWER option info for Gurobi.
%
%   DEFAULT_OPTS = MPOPTION_INFO_GUROBI('D')
%   VALID_OPTS   = MPOPTION_INFO_GUROBI('V')
%   EXCEPTIONS   = MPOPTION_INFO_GUROBI('E')
%
%   Returns a structure for Gurobi options for MATPOWER containing ...
%   (1) default options,
%   (2) valid options, or
%   (3) NESTED_STRUCT_COPY exceptions for setting options
%   ... depending on the value of the input argument.
%
%   This function is used by MPOPTION to set default options, check validity
%   of option names or modify option setting/copying behavior for this
%   subset of optional MATPOWER options.
%
%   See also MPOPTION.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
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

if nargin < 1
    selector = 'D';
end
if have_fcn('gurobi')
    switch upper(selector)
        case {'D', 'V'}     %% default and valid options
            opt = struct(...
                'gurobi',       struct(...
                    'method',       -1, ...
                    'timelimit',    Inf, ...
                    'threads',      0, ...
                    'opts',         [], ...
                    'opt_fname',    '', ...
                    'opt',          0 ...
                ) ...
            );
        case 'E'            %% exceptions used by nested_struct_copy() for applying
            opt = struct(...
                'name',         { 'gurobi.opts' }, ...
                'check',        0 ...
                );
%                 'copy_mode',    { @gurobi_options } ...
        otherwise
            error('mpoption_info_gurobi: ''%s'' is not a valid input argument', selector);
    end
else
    opt = struct([]);       %% Gurobi is not available
end
