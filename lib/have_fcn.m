function TorF = have_fcn(tag)
%HAVE_FCN  Test for optional functionality.
%   TorF = have_fcn(tag) returns 1 if the optional functionality is
%   available, 0 otherwise.
%
%   Possible values for input tag and their meanings:
%       bpmpd       - bp(), BPMPD interior point solver
%       constr      - constr(), solver from Optimization Toolbox 1.x/2.x
%       fmincon     - fmincon(), solver from Optimization Toolbox 2.x +
%       linprog     - linprog(), LP solver from Optimization Toolbox 2.x +
%       lp          - lp(), LP solver from Optimization Toolbox 1.x/2.x
%       minopf      - minopf(), MINOPF, MINOS-based OPF solver
%       quadprog    - quadprog(), QP solver from Optimization Toolbox 2.x +
%       qp          - qp(), QP solver from Optimization Toolbox 1.x/2.x
%       pdipmopf    - pdipmopf(), primal-dual interior point method OPF solver
%       scpdipmopf  - scpdipmopf(), step-controlled PDIPM OPF solver
%       smartmarket - runmkt() and friends, for running an auction
%       sparse_lp   - LP solver used by mp_lp can handle sparse matrices
%       sparse_qp   - QP solver used by mp_qp can handle sparse matrices
%       tralmopf    - tralmopf(), trust region based augmented Langrangian
%                     OPF solver

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

switch tag
    case 'bpmpd'
        TorF = exist(which('bp')) == 3;
    case 'constr'
        TorF = exist(which('constr')) == 2 && exist('foptions');
    case 'fmincon'
        TorF = exist(which('fmincon')) == 2;
    case 'linprog'
        TorF = exist(which('linprog')) == 2;
    case 'lp'
        TorF = exist(which('lp')) == 2;
    case 'minopf'
        TorF = exist(which('minopf')) == 3;
    case 'quadprog'
        TorF = exist(which('quadprog')) == 2;
    case 'qp'
        TorF = exist(which('qp')) == 2;
    case 'smartmarket'
        TorF = exist(which('runmkt')) == 2;
    case 'sparse_lp'
        TorF = have_fcn('bpmpd') | have_fcn('linprog');
    case 'sparse_qp'
        TorF = have_fcn('bpmpd') | have_fcn('quadprog');
    case {'pdipmopf', 'scpdipmopf', 'tralmopf'}
        v = ver('Matlab');
        %% requires >= Matlab 6.5 (R13) (released 20-Jun-2002)
        %% older versions do not have mxCreateDoubleScalar() function
        %% (they have mxCreateScalarDouble() instead)
        if datenum(v.Date) >= 731387
            switch tag
                case 'pdipmopf'
                    TorF = exist(which('pdipmopf')) == 3;
                case 'scpdipmopf'
                    TorF = exist(which('scpdipmopf')) == 3;
                case 'tralmopf'
                    %% requires >= Matlab 7.3 (R2006b) (released 03-Aug-2006)
                    %% older versions do not include the needed form of chol()
                    if datenum(v.Date) >= 732892
                        TorF = exist(which('tralmopf')) == 3;
                    else
                        TorF = 0;
                    end
            end
        else
            TorF = 0;
        end
    otherwise
        error('have_fcn: unknown functionality %s', tag);
end

return;
