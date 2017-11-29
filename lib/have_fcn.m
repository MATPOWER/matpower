function rv = have_fcn(tag, rtype)
%HAVE_FCN  Test for optional functionality / version info.
%   TORF = HAVE_FCN(TAG)
%   TORF = HAVE_FCN(TAG, TOGGLE)
%   VER_STR = HAVE_FCN(TAG, 'vstr')
%   VER_NUM = HAVE_FCN(TAG, 'vnum')
%   DATE    = HAVE_FCN(TAG, 'date')
%   INFO    = HAVE_FCN(TAG, 'all')
%
%   Returns availability, version and release information for optional
%   MATPOWER functionality. All information is cached, and the cached values
%   returned on subsequent calls. If the functionality exists, an attempt is
%   made to determine the release date and version number. The second
%   argument defines which value is returned, as follows:
%       <none>      1 = optional functionality is available, 0 = not available
%       'vstr'      version number as a string (e.g. '3.11.4')
%       'vnum'      version number as numeric value (e.g. 3.011004)
%       'date'      release date as a string (e.g. '20-Jan-2015')
%       'all'       struct with fields named 'av' (for 'availability'), 'vstr',
%                   'vnum' and 'date', and values corresponding to the above,
%                   respectively.
%
%   For functionality that is not available, all calls with a string-valued
%   second argument will return an empty value.
%
%   Alternatively, the optional functionality specified by TAG can be toggled
%   OFF or ON by calling HAVE_FCN with a numeric second argument TOGGLE with
%   one of the following values:
%       0 - turn OFF the optional functionality
%       1 - turn ON the optional functionality (if available)
%      -1 - toggle the ON/OFF state of the optional functionality
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point solver
%       clp         - CLP, LP/QP solver(http://www.coin-or.org/projects/Clp.xml)
%        opti_clp   -   version of CLP distributed with OPTI Toolbox
%                       (http://www.i2c2.aut.ac.nz/Wiki/OPTI/)
%       cplex       - CPLEX, IBM ILOG CPLEX Optimizer
%       e4st        - E4ST (http://e4st.com/)
%       fmincon     - FMINCON, solver from Optimization Toolbox 2.x +
%       fmincon_ipm - FMINCON with Interior Point solver, from Opt Tbx 4.x +
%       glpk        - GLPK, GNU Linear Programming Kit
%       gurobi      - GUROBI, Gurobi solver (http://www.gurobi.com/), 5.x +
%       intlinprog  - INTLINPROG, MILP solver from Optimization
%                     Toolbox 7.0 (R2014a)+
%       ipopt       - IPOPT, NLP solver
%                       (http://www.coin-or.org/projects/Ipopt.xml)
%       linprog     - LINPROG, LP solver from Optimization Toolbox 2.x +
%       linprog_ds  - LINPROG with dual-simplex solver
%                       from Optimization Toolbox 7.1 (R2014b) +
%       knitro      - KNITRO, NLP solver (http://www.ziena.com/)
%         knitromatlab - KNITRO, version 9.0.0+
%         ktrlink      - KNITRO, version < 9.0.0 (requires Opt Tbx)
%       matlab      - code is running under MATLAB, as opposed to Octave
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       most        - MOST, MATPOWER Optimal Scheduling Tool
%       mosek       - MOSEK, LP/QP solver (http://www.mosek.com/)
%       optimoptions - OPTIMOPTIONS, option setting funciton for Optim Tbx 6.3+
%       pardiso     - PARDISO, Parallel Sparse Direct & Iterative Linear Solver
%                       (http://www.pardiso-project.org)
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%       quadprog_ls - QUADPROG with large-scale interior point convex solver
%                       from Optimization Toolbox 6.x +
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       smartmarket - RUNMARKET and friends, for running an auction
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%       octave      - code is running under Octave, as opposed to MATLAB
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
%
%   Optional functionality can also be toggled OFF and ON by calling HAVE_FCN
%   with the following syntax,
%       TORF = HAVE_FCN(TAG, TOGGLE)
%   where TOGGLE takes a numeric value as follows:

%   Private tags for internal use only:
%       catchme         - support for 'catch me' syntax in try/catch constructs
%       evalc           - support for evalc() function
%       ipopt_auxdata   - support for ipopt_auxdata(), required by 3.11 and later
%       lu_vec          - support for lu(..., 'vector') syntax
%       pardiso_legacy  - PARDISO v5, individual MEX files for factor, solve, etc
%       pardiso_object  - PARDISO v6 and later, object interface
%       regexp_split    - support for 'split' argument to regexp()

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin > 1 && isnumeric(rtype)
    toggle = 1;
    on_off = rtype;
    if on_off < 0
        TorF = have_fcn(tag);
        on_off = ~TorF;
    end
else
    toggle = 0;
end

persistent fcns;

if toggle   %% change availability
    if on_off       %% turn on if available
        fcns = rmfield(fcns, tag);  %% delete field to force re-check
    else            %% turn off
        if ~isfield(fcns, tag)      %% not yet been checked
            TorF = have_fcn(tag);   %% cache result first
        end
        fcns.(tag).av = 0;          %% then turn off
    end
    TorF = have_fcn(tag);           %% return cached value
else        %% detect availability
    %% info not yet cached?
    if ~isfield(fcns, tag)
        %%-----  determine installation status, version number, etc.  -----
        %% initialize default values
        TorF = 0;
        vstr = '';
        rdate = '';

        switch tag
            %%-----  public tags  -----
            case 'bpmpd'
                TorF = exist('bp', 'file') == 3;
                if TorF
                    v = bpver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'clp'
                tmp = have_fcn('opti_clp', 'all');
                if tmp.av   %% have opti_clp
                    TorF = tmp.av;
                    vstr = tmp.vstr;
                    rdate = tmp.date;
                elseif exist('clp','file') == 2 && exist('mexclp','file') == 3
                    TorF = 1;
                    vstr = '';
                end
            case 'opti_clp'
                TorF = exist('opti_clp', 'file') == 2 && exist('clp', 'file') == 3;
                if TorF
                    str = evalc('clp');
                    pat = 'CLP: COIN-OR Linear Programming \[v([^\s,]+), Built ([^\],])+(,[^\]]*)*\]';  %% OPTI, Giorgetti/Currie
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        vstr = t{1}{1};
                        rdate = datestr(t{1}{2}, 'dd-mmm-yyyy');
                    end
                end
            case 'cplex'
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
                if TorF
                    try
                        cplex = Cplex('null');
                        vstr = cplex.getVersion;
                    catch
                        TorF = 0;
                    end
                end
            case 'e4st'
                TorF = exist('e4st_ver', 'file') == 2;
                if TorF
                    v = e4st_ver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case {'fmincon', 'fmincon_ipm', 'intlinprog', 'linprog', ...
                        'linprog_ds', 'optimoptions', 'quadprog', 'quadprog_ls'}
                matlab = have_fcn('matlab');
                if ~matlab || (matlab && license('test', 'optimization_toolbox'))
                    v = ver('optim');
                    if length(v) > 1
                        warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''optim'' on your path. Check each element of the output of ver(''optim'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
                        v = v(1);
                    end
                    if isempty(v) || isempty(v.Version)
                        vstr = '';      %% make sure it's a string
                        rdate = [];
                    else
                        vstr = v.Version;
                        rdate = v.Date;
                    end
                    otver = vstr2num(vstr);
                    switch tag
                        case 'fmincon'
                            TorF = (exist('fmincon', 'file') == 2 || ...
                                exist('fmincon', 'file') == 6) & matlab;
                        case 'intlinprog'
                            TorF = exist('intlinprog', 'file') == 2 & matlab;
                        case 'linprog'
                            TorF = exist('linprog', 'file') == 2 & matlab;  %% don't try to use Octave linprog
                        case 'quadprog'
                            TorF = exist('quadprog', 'file') == 2;
                            %% Octave optim 1.5.0 and earlier, had problems with
                            %% incorrect lambdas, including opposite sign
                            %% convention for equality multipliers
                            if ~matlab && otver <= 1.005
                                TorF = 0;
                            end
                        otherwise
                            if matlab
                                switch tag
                                    case 'fmincon_ipm'
                                        if otver >= 4       %% Opt Tbx 4.0+ (R208a+, MATLAB 7.6+)
                                            TorF = 1;
                                        end
                                    case 'linprog_ds'
                                        if otver >= 7.001   %% Opt Tbx 7.1+ (R2014b+, MATLAB 8.4+)
                                            TorF = 1;
                                        end
                                    case 'optimoptions'
                                        if otver >= 6.003   %% Opt Tbx 6.3+ (R2013a+, MATLAB 8.1+)
                                            TorF = 1;
                                        end
                                    case 'quadprog_ls'
                                        if otver >= 6       %% Opt Tbx 6.0+ (R2011a+, MATLAB 7.12+)
                                            TorF = 1;
                                        end
                                end
                            else    %% octave
                                TorF = 0;
                            end
                    end
                end
            case 'glpk'
                if exist('glpk','file') == 3    %% Windows OPTI install (no glpk.m)
                    TorF = 1;
                    str = evalc('glpk');
                    pat = 'GLPK: GNU Linear Programming Kit \[v([^\s,\]]+).*\]';  %% OPTI, Giorgetti/Currie
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        vstr = t{1}{1};
                    end
                    pat = 'Built ([^\],])+';  %% OPTI, Giorgetti/Currie
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        rdate = datestr(t{1}{1}, 'dd-mmm-yyyy');
                    end
                elseif exist('glpk','file') == 2    %% others have glpk.m and ...
                    if exist('__glpk__','file') == 3    %% octave __glpk__ MEX
                        TorF = 1;
                        if have_fcn('evalc')
                            str = evalc('glpk(1, 1, 1, 1, 1, ''U'', ''C'', -1, struct(''msglev'', 3))');
                            pat = 'GLPK Simplex Optimizer, v([^\s,]+)';
                            [s,e,tE,m,t] = regexp(str, pat);
                            if ~isempty(t)
                                vstr = t{1}{1};
                            end
                        end
                    elseif exist('glpkcc','file') == 3  %% MATLAB glpkcc MEX
                        TorF = 1;
                        str = evalc('glpk');
                        pat = 'GLPK Matlab interface\. Version: ([^\s,]+)';     %% glpkccm, Giorgetti/Klitgord
                        [s,e,tE,m,t] = regexp(str, pat);
                        if ~isempty(t)
                            vstr = t{1}{1};
                        end
                    end
                end
            case 'gurobi'
                TorF = exist('gurobi', 'file') == 3;
                if TorF
                    try
                        model = struct( ...
                            'A', sparse(1), ...
                            'rhs', 1, ...
                            'sense', '=', ...
                            'vtype', 'C', ...
                            'obj', 1, ...
                            'modelsense', 'min' ...
                        );
                        params = struct( ...
                            'outputflag', 0 ...
                        );
                        result = gurobi(model, params);
                        vstr = sprintf('%d.%d.%d', result.versioninfo.major, result.versioninfo.minor, result.versioninfo.technical);
                    catch % gurobiError
                        fprintf('Gurobi Error!\n');
%                         disp(gurobiError.message);
                    end
                end
            case 'ipopt'
                TorF = exist('ipopt', 'file') == 3;
                if TorF
                    str = evalc('qps_ipopt([],[1; 1],[1 1],[2],[2],[1; 1],[1; 1],[1; 1],struct(''verbose'', 2))');
                    pat = 'Ipopt version ([^\s,]+)';
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        vstr = t{1}{1};
                        if vstr2num(vstr) >= 3.011 && ...
                                ~exist('ipopt_auxdata', 'file')
                            TorF = 0;
                            warning('Improper installation of IPOPT. Version %s detected, but IPOPT_AUXDATA.M is missing.', vstr);
                        end
                    end
                end
            case 'knitro'       %% any Knitro
                tmp = have_fcn('knitromatlab', 'all');
                if tmp.av
                    TorF = tmp.av;
                    vstr = tmp.vstr;
                    rdate = tmp.date;
                else
                    tmp = have_fcn('ktrlink', 'all');
                    if tmp.av
                        TorF = tmp.av;
                        vstr = tmp.vstr;
                        rdate = tmp.date;
                    end
                end
            case {'knitromatlab', 'ktrlink'}
                %% knitromatlab for Knitro 9.0 or greater
                %% ktrlink for pre-Knitro 9.0, requires Optim Toolbox
                TorF = exist(tag, 'file') == 2;
                if TorF
                    try
                        str = evalc(['[x fval] = ' tag '(@(x)1,1);']);
                    end
                    TorF = exist('fval', 'var') && fval == 1;
                    if TorF
                        pat = 'KNITRO ([^\s]+)\n|Knitro ([^\s]+)\n';
                        [s,e,tE,m,t] = regexp(str, pat);
                        if ~isempty(t)
                            vstr = t{1}{1};
                        end
                    end
                end
            case 'matlab'
                v = ver('matlab');
                if length(v) > 1
                    warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
                    v = v(1);
                end
                if ~isempty(v) && isfield(v, 'Version') && ~isempty(v.Version)
                    TorF = 1;
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'minopf'
                TorF = exist('minopf', 'file') == 3;
                if TorF
                    v = minopfver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'most'
                TorF = exist('most', 'file') == 2;
                if TorF
                    v = mostver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'mosek'
                TorF = exist('mosekopt', 'file') == 3;
                if TorF
                    % MOSEK Version 6.0.0.93 (Build date: 2010-10-26 13:03:27)
                    % MOSEK Version 6.0.0.106 (Build date: 2011-3-17 10:46:54)
                    % MOSEK Version 7.0.0.134 (Build date: 2014-10-2 11:10:02)
                    pat = 'Version (\.*\d)+.*Build date: (\d+-\d+-\d+)';
                    [s,e,tE,m,t] = regexp(evalc('mosekopt'), pat);
                    if isempty(t)
                        [r, res] = mosekopt('version');
                        v = res.version;
                        vstr = sprintf('%d.%d.%d.%d', ...
                            v.major, v.minor, v.build, v.revision);
                    else
                        vstr = t{1}{1};
                        rdate = datestr(t{1}{2}, 'dd-mmm-yyyy');
                    end
                end
            case 'smartmarket'
                TorF = exist('runmarket', 'file') == 2;
                if TorF
                    v = mpver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'octave'
                TorF = exist('OCTAVE_VERSION', 'builtin') == 5;
                if TorF
                    v = ver('octave');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'pardiso'
                TorF = have_fcn('pardiso_object') || have_fcn('pardiso_legacy');
            case {'pdipmopf', 'scpdipmopf', 'tralmopf'}
                if have_fcn('matlab')
                    vn = have_fcn('matlab', 'vnum');
                    %% requires >= MATLAB 6.5 (R13) (released 20-Jun-2002)
                    %% older versions do not have mxCreateDoubleScalar() function
                    %% (they have mxCreateScalarDouble() instead)
                    if vn >= 6.005
                        switch tag
                            case 'pdipmopf'
                                TorF = exist('pdipmopf', 'file') == 3;
                            case 'scpdipmopf'
                                TorF = exist('scpdipmopf', 'file') == 3;
                            case 'tralmopf'
                                %% requires >= MATLAB 7.3 (R2006b) (released 03-Aug-2006)
                                %% older versions do not include the needed form of chol()
                                if vn >= 7.003
                                    TorF = exist('tralmopf', 'file') == 3;
                                end
                        end
                    end
                    if TorF
                        v = feval([tag 'ver'], 'all');
                        vstr = v.Version;
                        rdate = v.Date;
                    end
                end
            case 'sdp_pf'
                TorF = have_fcn('yalmip') && exist('mpoption_info_sdp_pf', 'file') == 2;
                if TorF
                    v = sdp_pf_ver('all');
                    vstr = v.Version;
                    rdate = v.Date;
                end
            case 'yalmip'
                TorF = ~have_fcn('octave') && exist('yalmip','file') == 2;
                %% YALMIP does not yet work with Octave, rdz 1/6/14
                if TorF
                    vstr = yalmip('version');
                    if length(vstr) == 8
                        yr = str2num(vstr(1:4));
                        mo = str2num(vstr(5:6));
                        dy = str2num(vstr(7:8));
                        rdate = datestr([yr mo dy 0 0 0], 'dd-mmm-yyyy');
                    end
                end
            case 'sdpt3'
                TorF = exist('sdpt3','file') == 2;
                if TorF
                    str = evalc('help sdpt3');
                    pat = 'version\s+([^\s]+).*Last Modified: ([^\n]+)\n';
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        vstr = t{1}{1};
                        rdate = datestr(t{1}{2}, 'dd-mmm-yyyy');
                    end
                end
            case 'sedumi'
                TorF = exist('sedumi','file') == 2;
                if TorF
                    warn_state = warning;  %% sedumi turns (and leaves!) off all warnings
                    str = evalc('x = sedumi([1 1], 1, [1;2])');
                    warning(warn_state);
                    pat = 'SeDuMi\s+([^\s]+)';
                    [s,e,tE,m,t] = regexp(str, pat);
                    if ~isempty(t)
                        vstr = t{1}{1};
                    end
                end

            %%-----  private tags  -----
            case 'catchme'  %% not supported by MATLAB <= 7.4 (R2007a), Octave <= 3.6
                if have_fcn('octave')
                    if have_fcn('octave', 'vnum') > 3.006
                        TorF = 1;
                    end
                else
                    if have_fcn('matlab', 'vnum') > 7.004
                        TorF = 1;
                    end
                end
            case 'evalc'
                if have_fcn('matlab')
                    TorF = 1;
                end
            case 'ipopt_auxdata'
                if have_fcn('ipopt')
                    vn = have_fcn('ipopt', 'vnum');
                    if ~isempty(vn) && vn >= 3.011
                        TorF = 1;
                    end
                end
            case 'lu_vec'       %% lu(..., 'vector') syntax supported?
                if have_fcn('matlab') && have_fcn('matlab', 'vnum') < 7.003
                    TorF = 0;     %% lu(..., 'vector') syntax not supported
                else
                    TorF = 1;
                end
            case 'pardiso_object'
                TorF = exist('pardiso', 'file') == 2;
                if TorF
                    try
                        id = 1;
                        A = sparse([1 2; 3 4]);
                        b = [1;1];
                        p = pardiso(id, 11, 0);
                        p.factorize(id, A);
                        x = p.solve(id, A, b);
                        p.free(id);
                        p.clear();
                        if any(x ~= [-1; 1])
                            TorF = 0;
                        end
                    catch
                        TorF = 0;
                    end
                end
            case 'pardiso_legacy'
                TorF = exist('pardisoinit', 'file') == 3 && ...
                        exist('pardisoreorder', 'file') == 3 && ...
                        exist('pardisofactor', 'file') == 3 && ...
                        exist('pardisosolve', 'file') == 3 && ...
                        exist('pardisofree', 'file') == 3;
                if TorF
                    try
                        A = sparse([1 2; 3 4]);
                        b = [1;1];
                        info = pardisoinit(11, 0);
                        info = pardisoreorder(A, info, false);
%                         % Summary PARDISO 5.1.0: ( reorder to reorder )
%                         pat = 'Summary PARDISO (\.*\d)+:';
%                         [s,e,tE,m,t] = regexp(evalc('info = pardisoreorder(A, info, true);'), pat);
%                         if ~isempty(t)
%                             vstr = t{1}{1};
%                         end
                        info = pardisofactor(A, info, false);
                        [x, info] = pardisosolve(A, b, info, false);
                        pardisofree(info);
                        if any(x ~= [-1; 1])
                            TorF = 0;
                        end
                    catch
                        TorF = 0;
                    end
                end
            case 'regexp_split'     %% missing for MATLAB < 7.3 & Octave < 3.8
                if have_fcn('matlab') && have_fcn('matlab', 'vnum') >= 7.003
                    TorF = 1;
                elseif have_fcn('octave', 'vnum') >= 3.008
                    TorF = 1;
                end

        %%-----  unknown tag  -----
            otherwise
                warning('have_fcn: unknown functionality ''%s''', tag);
                vstr = 'unknown';
        end

        %% assign values to cache
        fcns.(tag).av   = TorF;
        fcns.(tag).vstr = vstr;
        if isempty(vstr)
            fcns.(tag).vnum = [];
        else
            fcns.(tag).vnum = vstr2num(vstr);
        end
        fcns.(tag).date = rdate;
    end
end

%% extract desired values from cache
if nargin < 2 || toggle
    rv = fcns.(tag).av;
else
    switch lower(rtype)
        case 'vstr'
            rv = fcns.(tag).vstr;
        case 'vnum'
            rv = fcns.(tag).vnum;
        case 'date'
            rv = fcns.(tag).date;
        case 'all'
            rv = fcns.(tag);
    end
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
