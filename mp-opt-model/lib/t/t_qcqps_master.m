function t_qcqps_master(quiet)
% t_qcqps_master - Tests of QCQP solvers via qcqps_master.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%           1        2        3         4         5         6         7
algs = {'DEFAULT', 'MIPS', 'IPOPT', 'FMINCON', 'GUROBI', 'KNITRO', 'KNITRO_NLP'};
names = {'DEFAULT', 'MIPS', 'IPOPT', 'fmincon', 'GUROBI', 'KNITRO', 'KNITRO_NLP'};
check = {[], [], 'ipopt', 'fmincon', 'gurobi', 'knitro', 'knitro'};
%               1 2 3 4 5 6 7
does_lp      = [1 1 1 0 1 1 0];
does_qp      = [1 1 1 0 1 1 0];
does_nonconv = [1 1 1 1 1 1 1];

nlp = 16;
nqp = 28;
nqcqp_convex = 27;
nqcqp_nonconvex = 9;
nerrors = 13;
n = nlp+nqp+nqcqp_convex+nqcqp_nonconvex+nerrors;

t_begin(n * length(algs), quiet);

diff_alg_warn_id = 'optim:linprog:WillRunDiffAlg';
if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    s1 = warning('query', diff_alg_warn_id);
    warning('off', diff_alg_warn_id);
end

for k = 1:length(algs)
    if ~isempty(check{k}) && ~have_feature(check{k})
        t_skip(n, sprintf('%s not installed', names{k}));
    else
        opt = struct('verbose', 0, 'alg', algs{k});
        mpopt = struct( ...
            'verbose', 0, ...
            'opf', struct( ...
                'violation', 1e-6 ), ...
            'gurobi', struct( ...
                 'method', -1, ...
                 'timelimit', Inf, ...
                 'threads', 0, ...
                 'opts', [], ...
                 'opt_fname', [], ...
                 'opt', 0), ...
            'knitro', struct( ...
                 'tol_x', 1e-10, ...
                 'tol_f', 1e-10, ...
                 'maxit', 0, ...
                 'opt', 0) ...
        );
        opt.mips_opt.comptol = 1e-8;
        if have_feature('gurobi')
            opt.grb_opt = gurobi_options([], mpopt);
            opt.grb_opt.BarQCPConvTol = 1e-8;
        end
        if have_feature('knitro')
            opt.knitro_opt = artelys_knitro_options([],  mpopt);
            opt.knitro_opt.ncvx_qcqp_init = 1;
        end

        if does_lp(k)
            t = sprintf('%s - 3-d LP : ', names{k});
            %% 1) Based on example from 'doc linprog'
            b = [-5; -4; -6];
            A = [1 -1  1;
                -3  -2  -4;
                3  2  0];
            l = [-Inf; -42; -Inf];
            u = [20; Inf; 30];
            xmin = [0; 0; 0];
            xmax = [Inf; Inf; Inf];
            x0 = [];
            [x, f, s, out, lam] = qcqps_master([], b, [], [], [], [], A, l, u, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 15; 3], 6, [t 'x']);
            t_is(f, -78, 6, [t 'f']);
            t_is(lam.mu_l, [0;1.5;0], 7, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0;0.5], 7, [t 'lam.mu_u']);
            t_is(lam.lower, [1;0;0], 7, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 7, [t 'lam.upper']);
            
            %% 2) Same previous passing linear constraints via A and B
            t = sprintf('%s - 3-d LP using A and B: ', names{k});
            [x2, f2, s2, out2, lam2] = qcqps_master([], b, [], A(end,:), l(end), u(end), A(1:2,:), l(1:2), u(1:2), xmin, xmax, x0, opt);
            t_is(s2, s, 12, [t 'success']);
            t_is(x2, x, 12, [t 'x']);
            t_is(f2, f, 12, [t 'f']);
            t_is(lam.mu_l, lam2.mu_l, 12, [t 'lam.mu_l']);
            t_is(lam.mu_u, lam2.mu_u, 12, [t 'lam.mu_u']);
            t_is(lam.lower, lam2.lower, 12, [t 'lam.lower']);
            t_is(lam.upper, lam2.upper, 12, [t 'lam.upper']);
            t_ok(isequal(fieldnames(out),fieldnames(out2)), [t 'out']);            

            %% 3) Infeasible LP problem
            t = sprintf('%s - infeasible LP : ', names{k});
            p = struct('A', sparse([1 1]), 'b', [1;1], 'u', -1, 'xmin', [0;0], 'opt', opt);
            [x, f, s, out, lam] = qcqps_master(p);
            t_ok(s <= 0, [t 'no success']);
        else
            t_skip(nlp, sprintf('%s does not handle LP problems', names{k}));
        end

        if does_qp(k)
            t = sprintf('%s - unconstrained 3-d convex QP : ', names{k});
            %% 4) From http://www.akiti.ca/QuadProgEx0Constr.html
            H = [5 -2 -1; -2 4 3; -1 3 5];
            b = [2; -35; -47];
            x0 = [0; 0; 0];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], [], [], [], [], [], [], opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [3; 5; 7], 7, [t 'x']);
            t_is(f, -249, 11, [t 'f']);
            t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
            t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
            t_is(lam.lower, zeros(size(x)), 13, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

            t = sprintf('%s - constrained 2-d convex QP : ', names{k});
            %% 5) Example from 'doc quadprog'
            H = [   1   -1;
                    -1  2   ];
            b = [-2; -6];
            A = [   1   1;
                    -1  2;
                    2   1   ];
            l = [];
            u = [2; 2; 3];
            xmin = [0; 0];
            xmax = [Inf; Inf];
            x0 = [];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], A, l, u, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [2; 4]/3, 7, [t 'x']);
            t_is(f, -74/9, 6, [t 'f']);
            t_is(lam.mu_l, [0;0;0], 13, [t 'lam.mu_l']);
            t_is(lam.mu_u, [28;4;0]/9, 4, [t 'lam.mu_u']);
            t_is(lam.lower, zeros(size(x)), 7, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

            t = sprintf('%s - constrained 4-d convex QP : ', names{k});
            %% 6) From https://v8doc.sas.com/sashtml/iml/chap8/sect12.htm
            H = [   1003.1  4.3     6.3     5.9;
                    4.3     2.2     2.1     3.9;
                    6.3     2.1     3.5     4.8;
                    5.9     3.9     4.8     10  ];
            b = zeros(4,1);
            A = [   1       1       1       1;
                    0.17    0.11    0.10    0.18    ];
            l = [1; 0.10];
            u = [1; Inf];
            xmin = zeros(4,1);
            xmax = Inf(4,1);
            x0 = [1; 0; 0; 1];
            [x, f, s, out, lam] = qcqps_master(H, b, [], [], [], [], A, l, u, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);

            %% 7) Same previous passing a struct
            t = sprintf('%s - (struct) constrained 4-d convex QP : ', names{k});
            p = struct('H', H, 'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
            [x, f, s, out, lam] = qcqps_master(p);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 2.8; 0.2; 0]/3, 5, [t 'x']);
            t_is(f, 3.29/3, 6, [t 'f']);
            t_is(lam.mu_l, [6.58;0]/3, 6, [t 'lam.mu_l']);
            t_is(lam.mu_u, [0;0], 13, [t 'lam.mu_u']);
            t_is(lam.lower, [2.24;0;0;1.7667], 4, [t 'lam.lower']);
            t_is(lam.upper, zeros(size(x)), 13, [t 'lam.upper']);
        else
            t_skip(nqp, sprintf('%s does not handle QP problems', names{k}));
        end

        %% 8) From https://docs.gurobi.com/projects/examples/en/current/examples/matlab/qcp.html
        t = sprintf('%s - convex 3-d QCQP with lin objective: ', names{k});
        H = [];
        c = [-1;0;0];
        Q = cell(2,1);
        Q{1} = sparse([2 0 0; 0 2 0; 0 0 -2]);
        Q{2} = sparse([2 0 0; 0 0 -2; 0 -2 0]);
        B = zeros(2,3);
        lq = [-Inf;-Inf];
        uq = [0; 0];
        A = [1 1 1];
        l = 1;
        u = 1;
        xmin = zeros(3,1);
        xmax = Inf(3,1);
        x0 = zeros(3,1);
        [x, f, s, out, lam] = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [3.91577; 1.78203; 4.3022]*1e-1, 6, [t 'x']);
        t_is(f, -0.391577, 6, [t 'f']);
        t_is(lam.lower, [0.0482;0.0092;0.1658]*1e-8, 6, [t 'lam.lower']);
        t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
        t_is(lam.mu_l, 0, 6, [t 'lam.mu_l']);
        t_is(lam.mu_u, 0.391577, 6, [t 'lam.mu_u']);
        t_is(lam.mu_lq, [0; 0], 5, [t 'lam.mu_lq']);
        t_is(lam.mu_uq, [0.227544; 0.549342], 5, [t 'lam.mu_uq']);

        %% 9) Same previous passing a struct
        t = sprintf('%s - (struct) convex 3-d QCQP with lin objective : ', names{k});
        p = struct('H', H, 'c', c, 'Q', {Q}, 'B', B, 'lq', lq, 'uq', uq, ...
            'A', A, 'l', l, 'u', u, 'xmin', xmin, 'x0', x0, 'opt', opt);
        [x, f, s, out, lam] = qcqps_master(p);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [3.91577; 1.78203; 4.3022]*1e-1, 6, [t 'x']);
        t_is(f, -0.391577, 6, [t 'f']);
        t_is(lam.lower, [0.0482;0.0092;0.1658]*1e-8, 6, [t 'lam.lower']);
        t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
        t_is(lam.mu_l, 0, 6, [t 'lam.mu_l']);
        t_is(lam.mu_u, 0.391577, 6, [t 'lam.mu_u']);
        t_is(lam.mu_lq, [0; 0], 5, [t 'lam.mu_lq']);
        t_is(lam.mu_uq, [0.227544; 0.549342], 5, [t 'lam.mu_uq']);

        %% 10) From https://docs.mosek.com/latest/toolbox/examples-list.html#doc-example-file-qcqo1-m
        t = sprintf('%s - convex 3-d QCQP with quad objective: ', names{k});
        H = sparse([2 0 -1; 0 0.2 0; -1 0 2]);
        c = [0;-1;0];
        Q = {sparse([-2 0 0.2; 0 -2 0; 0.2 0 -0.2])};
        B = [1 1 1];
        lq = 1;
        uq = Inf;
        xmin = zeros(3,1);
        xmax = Inf(3,1);
        x0 = zeros(3,1);
        [x, f, s, out, lam] = qcqps_master(H, c, Q, B, lq, uq, [], [], [], xmin, xmax, x0, opt);
        t_is(s, 1, 12, [t 'success']);
        t_is(x, [4.4880; 9.3192; 6.7411]*1e-1, 4, [t 'x']);
        t_is(f, -0.4918, 4, [t 'f']);
        t_is(lam.lower, [0.2228;0.1073;0.1483]*1e-10, 6, [t 'lam.lower']);
        t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
        t_ok(isempty(lam.mu_l), [t 'lam.mu_l']);
        t_ok(isempty(lam.mu_u), [t 'lam.mu_u']);
        t_is(lam.mu_lq, 0.9419, 4, [t 'lam.mu_lq']);
        t_is(lam.mu_uq, 0, 4, [t 'lam.mu_uq']);

        %% 11) From "examples" folder of Knitro (exampleQCQP1)
        t = sprintf('%s - nonconvex 3-d QCQP : ', names{k});
        if does_nonconv(k)
            H = sparse([-2 -1 -1; -1 -4 0; -1 0 -2]);
            c = zeros(3,1);
            Q = {-2*speye(3)};
            B = sparse(1,3);
            lq = -Inf;
            uq = -25;
            A = [8 14 7];
            l = 56;
            u = 56;
            xmin = zeros(3,1);
            xmax = Inf(3,1);
            x0 = [0;0;20];
            [x, f, s, out, lam] = qcqps_master(H, c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [0; 0; 8], 5, [t 'x']);
            t_is(f, -64, 4, [t 'f']);
            if strcmp(algs{k}, 'GUROBI')
                %% See https://docs.gurobi.com/projects/optimizer/en/current/reference/attributes/constraintquadratic.html#qcpi
                t_skip(6, [t 'lam : GUROBI version 12.0 and earlier does not return multipliers for nonconvex QCQP']);
            else
                t_is(lam.lower, [10.285714;32;0], 5, [t 'lam.lower']);
                t_is(lam.upper, [0;0;0], 6, [t 'lam.upper']);
                t_is(lam.mu_l, 0, 5, [t 'lam.mu_l']);
                t_is(lam.mu_u, 2.28571, 5, [t 'lam.mu_u']);
                t_is(lam.mu_lq, 0, 4, [t 'lam.mu_lq']);
                t_is(lam.mu_uq, 0, 4, [t 'lam.mu_uq']);
            end
        else
            t_skip(nqcqp_nonconvex, sprintf('%s does not handle nonconvex QCQP problems', names{k}));
        end

        %% 12) Errors handling
        H = [];
        c = [-1;0;0];
        Q = cell(2,1);
        Q{1} = sparse([2 0 0; 0 2 0; 0 0 -2]);
        Q{2} = sparse([2 0 0; 0 0 -2; 0 -2 0]);
        H2 = Q{1};
        Q2 = Q; Q2{1} = Q{1}(2:end,2:end);
        B = zeros(2,3);
        B2 = zeros(3);
        lq = [-Inf;-Inf];
        uq = [0; 0];
        A = [1 1 1];
        l = 1;
        u = 1;
        xmin = zeros(3,1);
        xmax = Inf(3,1);
        x0 = zeros(3,1);

        t = sprintf('%s - error : ', names{k});
        try
            [x, f, s, out, lam] = qcqps_master([], c, [Q Q2], [], lq, uq, [], l, u, [], [], x0, opt);
        catch me
            msg = 'Q must be column vector cell array.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H2(:,1:2), c, Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'H must be a square matrix.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, Q2, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'All matrices Q{i}, i=1,...,2 must be square of the same size.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master([], [], Q, [], lq, uq, [], l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'Problem is incomplete: H, c, B and A can not be all empty.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H2, ones(1,2), Q, B, lq, uq, A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'Dimension of c (2) must be iqual to the number of variables (3).';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, Q, B2, lq, uq, A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'Dimension of B (3x3) should be number of quad constraints times number of variables (2x3).';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, [], B(:,2:end), lq, uq, A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'Dimension of B (2x2) should be number of quad constraints times number of variables (2x3).';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master([], c, [], [], lq, uq, [], l, u, [], [], x0, opt);
        catch me
            msg = 'No quadratic constraints were found. Inputs lq and uq should be empty.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, [], B, lq, uq, A(:,2:end), l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'The number of columns of matrix A (2) must be equal to the number of variables (3)';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, [], B, lq(2:end), uq(2:end), A, l, u, xmin, xmax, x0, opt);
        catch me
            msg = 'Dimension mismatch between uq, Q, and B.';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H, c, Q, B, lq, uq, A, lq, uq, xmin, xmax, x0, opt);
        catch me
            msg = 'Dimension of u (2) must be iqual to the number of linear constraints (1).';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H2, c, Q, B, lq, uq, A, l, u, xmin(2:end), xmax, x0, opt);
        catch me
            msg = 'Dimension of xmin (2) must be iqual to the number of variables (3).';
            t_ok(strfind(me.message, msg), [t msg]);
        end
        try
            [x, f, s, out, lam] = qcqps_master(H2, c, Q, B, lq, uq, A, l, u, xmin, xmax(2:end), x0, opt);
        catch me
            msg = 'Dimension of xmax (2) must be iqual to the number of variables (3).';
            t_ok(strfind(me.message, msg), [t msg]);
        end       
    end
end

if have_feature('quadprog') && have_feature('quadprog', 'vnum') == 7.005
    warning(s1.state, diff_alg_warn_id);
end

t_end;
