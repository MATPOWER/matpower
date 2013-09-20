function rv = mpver(varargin)
%MPVER  Prints or returns MATPOWER version info for current installation.
%   V = MPVER returns the current MATPOWER version number.
%   V = MPVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MPVER without assigning the
%   return value prints the version and release date of the current
%   installation of MATPOWER, MATLAB (or Octave), the Optimization Toolbox,
%   MIPS and any optional MATPOWER packages such as BPMPD_MEX, CPLEX,
%   Gurobi, Ipopt, Knitro, MINOPF, Mosek, PDIPMOPF, SCPDIPMOPF and TRAMLOPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005-2011 by Power System Engineering Research Center (PSERC)
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

%% the following only works if MATPOWER is explicitly in the path,
%% but not if it is only in the current working directory
% fs = filesep;
% p = fileparts(which('runpf'));
% if ~strcmp(p(1),fs)
%   [t, p] = strtok(p, filesep);
% end
% p = p(2:end);
% v{1} = ver(p);

v{1} = struct(  'Name',     'MATPOWER', ... 
                'Version',  '4.1+', ...
                'Release',  '', ...
                'Date',     '20-Sep-2013' );
if nargout > 0
    if nargin > 0
        rv = v{1};
    else
        rv = v{1}.Version;
    end
else
    if have_fcn('octave')
        v{2} = ver('octave');
    else
        v{2} = ver('matlab');
    end
    v{3} = ver('optim');
    for n = 1:3
        if n == 3 && isempty(v{3})
            fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
            continue;
        elseif n == 3 && ~license('test', 'optimization_toolbox')
            fprintf('\n%-22s -- no license --', 'Optimization Toolbox');
            continue;
        end
        fprintf('\n%-22s Version %-9s', v{n}.Name, v{n}.Version);
        if ~isempty(v{n}.Date)
            fprintf('  %11s', v{n}.Date);
            if ~isempty(v{n}.Release)
                fprintf('   Release: %-10s', v{n}.Release);
            end
        end
    end
    fprintf('\n');
    mipsver;
    if have_fcn('bpmpd')
        if exist('bpver', 'file') == 2
            bpver;
        else
            fprintf('%-22s Version 2.21 or earlier\n', 'BPMPD_MEX');
        end
    else
        fprintf('%-22s -- not installed --\n', 'BPMPD_MEX');
    end
    if have_fcn('cplex')
        cplex = Cplex('null');
        fprintf('%-22s Version %-10s %-11s   %s\n', 'CPLEX', cplex.getVersion, '', computer);
    else
        fprintf('%-22s -- not installed --\n', 'CPLEX');
    end
    if have_fcn('gurobi')
        gurobiver;
    else
        fprintf('%-22s -- not installed --\n', 'Gurobi');
    end
    if have_fcn('ipopt')
        str = evalc('qps_ipopt([],1,1,1,1,1,1,1,struct(''verbose'', 2))');
        pat = 'Ipopt version ([^\s,]+)';
        [s,e,tE,m,t] = regexp(str, pat);
        if isempty(t)
            vn = '<unknown>';
        else
            vn = t{1}{1};
        end
        fprintf('%-22s Version %-10s %-11s   %s\n', 'IPOPT', vn, '', computer);
    else
        fprintf('%-22s -- not installed --\n', 'IPOPT');
    end
    if have_fcn('knitro')
        str = evalc('[x fval] = ktrlink(@(x)1,1);');
        pat = 'KNITRO ([^\s]+)\n';
        [s,e,tE,m,t] = regexp(str, pat);
        if isempty(t)
            vn = '<unknown>';
        else
            vn = t{1}{1};
        end
        fprintf('%-22s Version %-10s %-11s   %s\n', 'KNITRO', vn, '', computer);
    else
        fprintf('%-22s -- not installed --\n', 'KNITRO');
    end
    if have_fcn('minopf')
        if exist('minopfver', 'file') == 2
            minopfver;
        else
            fprintf('%-22s Version 3.0b2 or earlier\n', 'MINOPF');
        end
    else
        fprintf('%-22s -- not installed --\n', 'MINOPF');
    end
    if have_fcn('mosek')
        % (this code is also in qps_mosek.m)
        % MOSEK Version 6.0.0.93 (Build date: 2010-10-26 13:03:27)
        % MOSEK Version 6.0.0.106 (Build date: 2011-3-17 10:46:54)
%        pat = 'Version (\.*\d)+.*Build date: (\d\d\d\d-\d\d-\d\d)';
        pat = 'Version (\.*\d)+.*Build date: (\d+-\d+-\d+)';
        [s,e,tE,m,t] = regexp(evalc('mosekopt'), pat);
        if isempty(t)
            vn = '<unknown>';
            d  = '';
        else
            vn = t{1}{1};
            d  = datestr(t{1}{2}, 'dd-mmm-yyyy');
        end
        fprintf('%-22s Version %-10s %-11s   %s\n', 'MOSEK', vn, d, computer);
    else
        fprintf('%-22s -- not installed --\n', 'MOSEK');
    end
    if have_fcn('pdipmopf')
        pdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'PDIPMOPF');
    end
    if have_fcn('scpdipmopf')
        scpdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'SCPDIPMOPF');
    end
    if have_fcn('tralmopf')
        tralmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'TRALMOPF');
    end

    fprintf('%-22s %s\n\n', 'Architecture:', computer);
    
    fprintf('  MATPOWER %s is distributed under the GNU General Public License.\n', v{1}.Version);
    fprintf('  Please see the LICENSE and COPYING files for details.\n\n');
end
