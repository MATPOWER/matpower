function TorF = have_fcn(tag)
%HAVE_FCN  Test for optional functionality.
%   TORF = HAVE_FCN(TAG) returns 1 if the optional functionality is
%   available, 0 otherwise.
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point solver
%       constr      - CONSTR, solver from Optimization Toolbox 1.x/2.x
%       cplex       - CPLEX, IBM ILOG CPLEX Optimizer
%       fmincon     - FMINCON, solver from Optimization Toolbox 2.x +
%       ipopt       - IPOPT, NLP solver (https://projects.coin-or.org/Ipopt/)
%       linprog     - LINPROG, LP solver from Optimization Toolbox 2.x +
%       lp          - LP, LP solver from Optimization Toolbox 1.x/2.x
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       mosek       - MOSEK, LP/QP solver (http://www.mosek.com/)
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%       qp          - QP, QP solver from Optimization Toolbox 1.x/2.x
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       smartmarket - RUNMARKET and friends, for running an auction
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%       anon_fcns   - anonymous functions, MATLAB version >= 7
%       octave      - code is running under Octave, not MATLAB
%
%   Examples:
%       if have_fcn('minopf')
%           results = runopf(mpc, mpoption('OPF_ALG', 500));
%       end

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
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

switch tag
    case 'bpmpd'
        TorF = exist('bp', 'file') == 3;
    case 'constr'
        TorF = exist('constr', 'file') == 2 && exist('foptions', 'file');
    case 'cplex'
        TorF = 0;
        if exist('cplexqp', 'file')
            %% it's installed, but we need to check for MEX for this arch
            p = which('cplexqp');   %% get the path
            len = length(p) - length('cplexqp.p');
            w = what(p(1:len));             %% look for mex files on the path
            for k = 1:length(w.mex)
                if regexp(w.mex{k}, 'cplexlink[^\.]*');
                    TorF = 1;
                    break;
                end
            end
        end
    case 'fmincon'
        TorF = exist('fmincon', 'file') == 2;
    case 'ipopt'
        TorF = exist('ipopt', 'file') == 3;
    case 'linprog'
        TorF = exist('linprog', 'file') == 2;
    case 'lp'
        TorF = exist('lp', 'file') == 2;
    case 'minopf'
        TorF = exist('minopf', 'file') == 3;
    case 'mosek'
        TorF = exist('mosekopt', 'file') == 3;
    case 'quadprog'
        TorF = exist('quadprog', 'file') == 2;
    case 'qp'
        TorF = exist('qp', 'file') == 2;
    case 'smartmarket'
        TorF = exist('runmarket', 'file') == 2;
    case 'octave'
        TorF = exist('OCTAVE_VERSION', 'builtin') == 5;
    case 'anon_fcns'
        if have_fcn('octave')
            TorF = 1;
        else
            v = ver('Matlab');
            if str2double(v.Version(1)) < 7    %% anonymous functions not available
                TorF = 0;
            else
                TorF = 1;
            end
        end
    case {'pdipmopf', 'scpdipmopf', 'tralmopf'}
        if have_fcn('octave')
            TorF = 0;
        else
            v = ver('Matlab');
            %% requires >= MATLAB 6.5 (R13) (released 20-Jun-2002)
            %% older versions do not have mxCreateDoubleScalar() function
            %% (they have mxCreateScalarDouble() instead)
            if datenum(v.Date) >= 731387
                switch tag
                    case 'pdipmopf'
                        TorF = exist('pdipmopf', 'file') == 3;
                    case 'scpdipmopf'
                        TorF = exist('scpdipmopf', 'file') == 3;
                    case 'tralmopf'
                        %% requires >= MATLAB 7.3 (R2006b) (released 03-Aug-2006)
                        %% older versions do not include the needed form of chol()
                        if datenum(v.Date) >= 732892
                            TorF = exist('tralmopf', 'file') == 3;
                        else
                            TorF = 0;
                        end
                end
            else
                TorF = 0;
            end
        end
    otherwise
        error('have_fcn: unknown functionality %s', tag);
end
