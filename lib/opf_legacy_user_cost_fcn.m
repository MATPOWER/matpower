function [f, df, d2f] = opf_legacy_user_cost_fcn(x, cp)
%OPF_LEGACY_USER_COST_FCN  Evaluates legacy user costs and derivatives.
%   [F, DF, D2F] = OPF_LEGACY_USER_COST_FCN(X, CP)
%
%   Evaluates the legacy user-defined costs and derivatives.
%
%   Inputs:
%     X : cell array with vectors of optimization variables
%     CP : legacy user-defined cost parameter struct such as returned by
%          @OPT_MODEL/GET_COST_PARAMS
%
%   Outputs:
%     F  : sum of generator polynomial costs
%     DF : (optional) gradient (column vector) of polynomial costs
%     D2F : (optional) Hessian of polynomial costs
%
%   Examples:
%       f = opf_legacy_user_cost_fcn(x, cp);
%       [f, df] = opf_legacy_user_cost_fcn(x, cp);
%       [f, df, d2f] = opf_legacy_user_cost_fcn(x, cp);

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% unpack data
[N, Cw, H, dd, rh, kk, mm] = ...
    deal(cp.N, cp.Cw, cp.H, cp.dd, cp.rh, cp.kk, cp.mm);
nx = length(x);

if isempty(N)
    f = 0;
    df = zeros(nx, 1);
    d2f = sparse(nx, nx);
else
    nw = size(N, 1);
    r = N * x - rh;                 %% Nx - rhat
    iLT = find(r < -kk);            %% below dead zone
    iEQ = find(r == 0 & kk == 0);   %% dead zone doesn't exist
    iGT = find(r > kk);             %% above dead zone
    iND = [iLT; iEQ; iGT];          %% rows that are Not in the Dead region
    iL = find(dd == 1);             %% rows using linear function
    iQ = find(dd == 2);             %% rows using quadratic function
    LL = sparse(iL, iL, 1, nw, nw);
    QQ = sparse(iQ, iQ, 1, nw, nw);
    kbar = sparse(iND, iND, [   ones(length(iLT), 1);
                                zeros(length(iEQ), 1);
                                -ones(length(iGT), 1)], nw, nw) * kk;
    rr = r + kbar;                  %% apply non-dead zone shift
    M = sparse(iND, iND, mm(iND), nw, nw);  %% dead zone or scale
    diagrr = sparse(1:nw, 1:nw, rr, nw, nw);
    
    %% linear rows multiplied by rr(i), quadratic rows by rr(i)^2
    w = M * (LL + QQ * diagrr) * rr;

    f = (w' * H * w) / 2 + Cw' * w;

    %%----- evaluate cost gradient -----
    if nargout > 1
        HwC = H * w + Cw;
        AA = N' * M * (LL + 2 * QQ * diagrr);
        df = AA * HwC;

        %% numerical check
        if 0    %% 1 to check, 0 to skip check
            ddff = zeros(size(df));
            step = 1e-7;
            tol  = 1e-3;
            for k = 1:length(x)
                xx = x;
                xx(k) = xx(k) + step;
                ddff(k) = (opf_costfcn(xx, om) - f) / step;
            end
            if max(abs(ddff - df)) > tol
                idx = find(abs(ddff - df) == max(abs(ddff - df)));
                fprintf('\nMismatch in gradient\n');
                fprintf('idx             df(num)         df              diff\n');
                fprintf('%4d%16g%16g%16g\n', [ 1:length(df); ddff'; df'; abs(ddff - df)' ]);
                fprintf('MAX\n');
                fprintf('%4d%16g%16g%16g\n', [ idx'; ddff(idx)'; df(idx)'; abs(ddff(idx) - df(idx))' ]);
                fprintf('\n');
            end
        end     %% numeric check

        %% ---- evaluate cost Hessian -----
        if nargout > 2
            d2f = AA * H * AA' + 2 * N' * M * QQ * sparse(1:nw, 1:nw, HwC, nw, nw) * N;
        end
    end
end
