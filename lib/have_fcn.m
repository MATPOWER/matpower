function TorF = have_fcn(tag)
%HAVE_FCN  Test for optional functionality.
%   TORF = HAVE_FCN(TAG) returns 1 if the optional functionality is
%   available, 0 otherwise.
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point solver
%       cplex       - CPLEX, IBM ILOG CPLEX Optimizer
%       fmincon     - FMINCON, solver from Optimization Toolbox 2.x +
%       fmincon_ipm - FMINCON with Interior Point solver, from Opt Tbx 4.x +
%       glpk        - GLPK, GNU Linear Programming Kit
%       gurobi      - GUROBI, Gurobi solver (http://www.gurobi.com/), 5.x +
%       ipopt       - IPOPT, NLP solver (https://projects.coin-or.org/Ipopt/)
%       linprog     - LINPROG, LP solver from Optimization Toolbox 2.x +
%       linprog_ds  - LINPROG with dual-simplex solver
%                       from Optimization Toolbox 7.1 (R2014b) +
%       knitro      - KNITRO, NLP solver (http://www.ziena.com/)
%         knitromatlab - KNITRO, version 9.0.0+
%         ktrlink      - KNITRO, version < 9.0.0 (requires Opt Tbx)
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       mosek       - MOSEK, LP/QP solver (http://www.mosek.com/)
%       optimoptions - OPTIMOPTIONS, option setting funciton for Optim Tbx 6.3+
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%       quadprog_ls - QUADPROG with large-scale interior point convex solver
%                       from Optimization Toolbox 6.x +
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       smartmarket - RUNMARKET and friends, for running an auction
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%       octave      - code is running under Octave, not MATLAB
%       sdp_pf      - SDP_PF applications of semi-definite programming
%                     relaxation of power flow equations
%       yalmip      - YALMIP SDP modeling platform
%       sedumi      - SeDuMi SDP solver
%       sdpt3       - SDPT3 SDP solver
%
%   Examples:
%       if have_fcn('minopf')
%           results = runopf(mpc, mpoption('opf.ac.solver', 'MINOPF'));
%       end

%   Private tags for internal use only:
%       catchme         - support for 'catch me' syntax in try/catch constructs
%       ipopt_auxdata   - support for ipopt_auxdata(), required by 3.11 and later
%       regexp_split    - support for 'split' argument to regexp()

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2014 by Power System Engineering Research Center (PSERC)
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

persistent fcns;

if isfield(fcns, tag) && ~isempty(fcns.(tag))
    TorF = fcns.(tag);
else
    switch tag
    %%-----  public tags  -----
        case 'bpmpd'
            TorF = exist('bp', 'file') == 3;
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
        case {'fmincon', 'fmincon_ipm', 'linprog', 'linprog_ds', ...
                    'optimoptions', 'quadprog', 'quadprog_ls'}
            if license('test', 'optimization_toolbox')
                switch tag
                    case 'fmincon'
                        TorF = exist('fmincon', 'file') == 2 || ...
                            exist('fmincon', 'file') == 6;
                    case 'linprog'
                        TorF = exist('linprog', 'file') == 2;
                    case 'quadprog'
                        TorF = exist('quadprog', 'file') == 2;
                    otherwise
                        v = ver('optim');
                        otver = vstr2num(v.Version);
                        switch tag
                            case 'fmincon_ipm'
                                if otver >= 4
                                    TorF = 1;
                                else
                                    TorF = 0;
                                end
                            case 'linprog_ds'
                                if otver >= 7.001
                                    TorF = 1;
                                else
                                    TorF = 0;
                                end
                            case 'optimoptions'
                                if otver >= 6.003
                                    TorF = 1;
                                else
                                    TorF = 0;
                                end
                            case 'quadprog_ls'
                                if otver >= 6
                                    TorF = 1;
                                else
                                    TorF = 0;
                                end
                        end
                end
            else
                TorF = 0;
            end
        case 'glpk'     %% Octave version must be later than 3.4, so check for 'catchme'
            TorF = have_fcn('catchme') && exist('glpk','file') == 2 && ...
                (exist('__glpk__','file') == 3 || exist('glpkcc','file') == 3);
        case 'gurobi'
            TorF = exist('gurobi', 'file') == 3;
        case 'ipopt'
            TorF = exist('ipopt', 'file') == 3;
        case 'knitro'       %% any Knitro
            TorF = have_fcn('knitromatlab') || have_fcn('ktrlink');
        case 'knitromatlab'     %% Knitro 9.0 or greater
            TorF = exist('knitromatlab', 'file') == 2;
            if TorF
                try
                    str = evalc('[x fval] = knitromatlab(@(x)1,1);');
                end
                TorF = exist('fval', 'var') && fval == 1;
            end
        case 'ktrlink'      %% pre 9.0 Knitro interface, requires Opt Tbx
            TorF = exist('ktrlink', 'file') == 2;
            if TorF
                try
                    str = evalc('[x fval] = ktrlink(@(x)1,1);');
                end
                TorF = exist('fval', 'var') && fval == 1;
            end
        case 'minopf'
            TorF = exist('minopf', 'file') == 3;
        case 'mosek'
            TorF = exist('mosekopt', 'file') == 3;
        case 'smartmarket'
            TorF = exist('runmarket', 'file') == 2;
        case 'octave'
            TorF = exist('OCTAVE_VERSION', 'builtin') == 5;
        case {'pdipmopf', 'scpdipmopf', 'tralmopf'}
            if have_fcn('octave')
                TorF = 0;
            else
                v = ver('Matlab');
                %% requires >= MATLAB 6.5 (R13) (released 20-Jun-2002)
                %% older versions do not have mxCreateDoubleScalar() function
                %% (they have mxCreateScalarDouble() instead)
                if vstr2num(v.Version) >= 6.005
                    switch tag
                        case 'pdipmopf'
                            TorF = exist('pdipmopf', 'file') == 3;
                        case 'scpdipmopf'
                            TorF = exist('scpdipmopf', 'file') == 3;
                        case 'tralmopf'
                            %% requires >= MATLAB 7.3 (R2006b) (released 03-Aug-2006)
                            %% older versions do not include the needed form of chol()
                            if vstr2num(v.Version) >= 7.003
                                TorF = exist('tralmopf', 'file') == 3;
                            else
                                TorF = 0;
                            end
                    end
                else
                    TorF = 0;
                end
            end
        case 'sdp_pf'
            TorF = have_fcn('yalmip') && exist('mpoption_info_sdp_pf', 'file') == 2;
        case 'yalmip'
            TorF = ~have_fcn('octave') && exist('yalmip','file') == 2;
            %% YALMIP does not yet work with Octave, rdz 1/6/14
        case 'sedumi'
            TorF = exist('sedumi','file') == 2;
        case 'sdpt3'
            TorF = exist('sdpt3','file') == 2;

    %%-----  private tags  -----
        case 'catchme'  %% not supported by Matlab <= 7.4 (R2007a), Octave <= 3.6
            if have_fcn('octave')
                v = ver('Octave');
                if vstr2num(v.Version) <= 3.006
                    TorF = 0;
                else
                    TorF = 1;
                end
            else
                v = ver('Matlab');
                if vstr2num(v.Version) <= 7.004
                    TorF = 0;
                else
                    TorF = 1;
                end
            end
        case 'ipopt_auxdata'
            if have_fcn('ipopt')
                str = evalc('qps_ipopt([],1,1,1,1,1,1,1,struct(''verbose'', 2))');
                pat = 'Ipopt version ([^\s,]+)';
                [s,e,tE,m,t] = regexp(str, pat);
                if isempty(t)
                    TorF = 0;       %% assume version is less than 3.11
                else
                    vn = t{1}{1};
                    if vstr2num(vn) >= 3.011
                        TorF = 1;
                    else
                        TorF = 0;
                    end
                end
            else
                TorF = 0;
            end
        case 'regexp_split'
            TorF = 1;
            if have_fcn('octave')   %% only missing for Octave < 3.8
                v = ver('Octave');
                if vstr2num(v.Version) < 3.008
                    TorF = 0;
                end
            end

    %%-----  unknown tag  -----
        otherwise
            error('have_fcn: unknown functionality %s', tag);
    end
    fcns.(tag) = TorF;
end


function num = vstr2num(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
