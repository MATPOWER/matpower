function TorF = have_fcn(tag)
%HAVE_FCN  Test for optional functionality.
%   TORF = HAVE_FCN(TAG) returns 1 if the optional functionality is
%   available, 0 otherwise.
%
%   Possible values for input TAG and their meanings:
%       bpmpd       - BP, BPMPD interior point solver
%       constr      - CONSTR, solver from Optimization Toolbox 1.x/2.x
%       fmincon     - FMINCON, solver from Optimization Toolbox 2.x +
%       linprog     - LINPROG, LP solver from Optimization Toolbox 2.x +
%       lp          - LP, LP solver from Optimization Toolbox 1.x/2.x
%       minopf      - MINOPF, MINOPF, MINOS-based OPF solver
%       quadprog    - QUADPROG, QP solver from Optimization Toolbox 2.x +
%       qp          - QP, QP solver from Optimization Toolbox 1.x/2.x
%       pdipmopf    - PDIPMOPF, primal-dual interior point method OPF solver
%       scpdipmopf  - SCPDIPMOPF, step-controlled PDIPM OPF solver
%       smartmarket - RUNMARKET and friends, for running an auction
%       tralmopf    - TRALMOPF, trust region based augmented Langrangian
%                     OPF solver
%       anon_fcns   - anonymous functions, Matlab version >= 7
%       octave      - code is running under Octave, not Matlab
%
%   Examples:
%       if have_fcn('minopf')
%           results = runopf(mpc, mpoption('OPF_ALG', 500));
%       end

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004-2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

switch tag
    case 'bpmpd'
        TorF = exist('bp', 'file') == 3;
    case 'constr'
        TorF = exist('constr', 'file') == 2 && exist('foptions', 'file') ...
            && ~have_fcn('octave');
    case 'fmincon'
        TorF = exist('fmincon', 'file') == 2;
    case 'linprog'
        TorF = exist('linprog', 'file') == 2;
    case 'lp'
        TorF = exist('lp', 'file') == 2;
    case 'minopf'
        TorF = exist('minopf', 'file') == 3;
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
            %% requires >= Matlab 6.5 (R13) (released 20-Jun-2002)
            %% older versions do not have mxCreateDoubleScalar() function
            %% (they have mxCreateScalarDouble() instead)
            if datenum(v.Date) >= 731387
                switch tag
                    case 'pdipmopf'
                        TorF = exist('pdipmopf', 'file') == 3;
                    case 'scpdipmopf'
                        TorF = exist('scpdipmopf', 'file') == 3;
                    case 'tralmopf'
                        %% requires >= Matlab 7.3 (R2006b) (released 03-Aug-2006)
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
