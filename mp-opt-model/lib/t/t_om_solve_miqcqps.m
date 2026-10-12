function t_om_solve_miqcqps(quiet)
% t_om_solve_miqcqps - Tests of MIQP/MIQCQP solvers via opt_model.solve.

%   MP-Opt-Model
%   Copyright (c) 2019-2026, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%           1         2         3        4
algs  = {'DEFAULT', 'CPLEX', 'GUROBI', 'MOSEK'};
names = {'DEFAULT', 'CPLEX', 'Gurobi', 'MOSEK'};
check = {@have_miqp_solver, 'cplex', 'gurobi', 'mosek'};

%               1                     2 3 4
does_miqp =    [have_miqp_solver()    1 1 1];
does_miqcqp =  [have_miqp_solver()    0 1 0];
does_nonconv = [have_feature('cplex') 0 0 0];
% replace w/lines below once CPLEX & MOSEK interfaces for MIQCQP are implemented
% does_miqcqp =  [have_miqp_solver()    1 1 1];
% does_nonconv = [have_feature('cplex') 1 0 0];

nmiqp_convex = 3;
nmiqp_nonconvex = 6; 
nmiqcqp_convex = 6;
nmiqcqp_nonconvex = 9; 

n = nmiqp_convex+nmiqp_nonconvex+nmiqcqp_convex+nmiqcqp_nonconvex;

t_begin(n * length(algs), quiet);

for k = 1:length(algs)
    if ~isempty(check{k}) && ...
            (ischar(check{k}) && ~have_feature(check{k}) || ...
             isa(check{k}, 'function_handle') && ~check{k}())
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
                 'opt', 0) ...
        );
        if have_feature('gurobi')
            opt.grb_opt = gurobi_options([], mpopt);
            opt.grb_opt.BarQCPConvTol = 1e-8;
            opt.gurobi.MIPGap = 1e-9;
            opt.gurobi.MIPGapAbs = 1e-9;
            opt.gurobi.FeasibilityTol = 1e-9;
            opt.gurobi.OptimalityTol = 1e-9;
            opt.gurobi.IntFeasTol  = 1e-9;
        end
        opt_r = opt;
        opt_r.relax_integer = 1;
        opt_f = opt;
        opt_f.fix_integer = 1;

        if does_miqp(k)
            t = sprintf('%s - convex 2-d MIQP : ', names{k});
            %% 1) From OPTI Toolbox: https://jonathancurrie.github.io/OPTI/examples/problem-types/miqp/
            H = [1  -1; -1  2];
            c = [-2 -6]';
            Q = {spalloc(2,2,0)}; B = []; lq = []; uq = [];
            A = [1  1; -1  2; 2 1];
            u = [2 2 3]';
            l = [];
            xmin = zeros(2,1);
            xmax = Inf(2,1);
            x0 = [];
            vtype = 'IC';
            om = opt_model();
            om.init_set_types();
            om.var.add('x', 2, x0, xmin, xmax, vtype);
            om.qdc.add(om.var, 'cost', H, c);
            om.qcn.add(om.var, 'xQx', Q, B, lq, uq);
            om.lin.add(om.var, 'Ax', A, l, u);
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [1, 1]', 7, [t 'x']);
            t_is(f, -7.500, 4, [t 'f']);
            
            %% 2) From Billionnet et al. 2012 (https://hal.science/hal-01125718v1/file/Extending_the_QCR_method_to_general_mixed-integer_.pdf)
            t = sprintf('%s - nonconvex 4-d MIQP : ', names{k});
            if does_nonconv(k)
                H = [ -14    6  -30   -8;
                        6  -28  -14  -26;
                      -30  -14   16   14;
                       -8  -26   14   24 ];
                c = [15; 10; -7; -4];
                Q = {spalloc(4,4,0)}; B  = []; lq = []; uq = [];
                A = [5 1 8 4];
                u = 95;
                l = [];
                xmin = zeros(4,1);
                xmax = 10*ones(4,1);
                x0 = [];
                vtype = 'IICC';
                om = opt_model();
                om.init_set_types();
                om.var.add('x', 4, x0, xmin, xmax, vtype);
                om.qdc.add(om.var, 'cost', H, c);
                om.qcn.add(om.var, 'xQx', Q, B, lq, uq);
                om.lin.add(om.var, 'Ax', A, l, u);
                [x, f, s, out, lam] = om.solve(opt);
                t_is(s, 1, 12, [t 'success']);
                t_is(x, [8, 10, 2.0268, 7.1964]', 2, [t 'x']);
                t_is(f, -3434.2700, 2, [t 'f']);

                %% 3) Same previous fixing integer variables to their optimal values
                t = sprintf('%s - nonconvex 4-d MIQP (integer fixed) : ', names{k});
                x0 = [8, 10, 0, 0]';
                om.var.set_params('x','v0',x0);
                [x, f, s, out, lam] = om.solve(opt_f);
                t_is(s, 1, 12, [t 'success']);
                t_is(x, [8, 10, 2.0268, 7.1964]', 2, [t 'x']);
                t_is(f, -3434.2700, 2, [t 'f']);
            else
                t_skip(nmiqp_nonconvex, sprintf('%s does not handle nonconvex MIQP problems', names{k}));
            end
        else
            t_skip(nmiqp_convex+nmiqp_nonconvex, sprintf('%s does not handle MIQP problems', names{k}));
        end
        
        %% 4) From OPTI Toolbox: https://jonathancurrie.github.io/OPTI/examples/problem-types/miqcqp/
        if does_miqcqp(k)
            t = sprintf('%s - convex 2-d MIQCQP : ', names{k});
            H = [1  0; 0  1];
            c = [-2 -2]';
            Q = {2*speye(2)};
            B  = [0 -2];
            lq = [];
            uq = 1;
            A = [-1 1; 1 3];
            u = [2 5]';
            l = [];
            xmin = zeros(2,1);
            xmax = Inf(2,1);
            x0 = [];
            vtype = 'IC';
            om = opt_model();
            om.init_set_types();
            om.var.add('x', 2, x0, xmin, xmax, vtype);
            om.qdc.add(om.var, 'cost', H, c);
            om.qcn.add(om.var, 'xQx', Q, B, lq, uq);
            om.lin.add(om.var, 'Ax', A, l, u);
            [x, f, s, out, lam] = om.solve(opt);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [1.0, 4/3]', 7, [t 'x']);
            t_is(f, -3.2778, 4, [t 'f']);

            %% 5) Same previous relaxing integer variable
            t = sprintf('%s - convex 2-d MIQCQP (integer relaxed) : ', names{k});
            [x, f, s, out, lam] = om.solve(opt_r);
            t_is(s, 1, 12, [t 'success']);
            t_is(x, [7/5, 6/5]', 7, [t 'x']);
            t_is(f, -3.5000, 4, [t 'f']);

            %% 6) From You and Dai 2020 (https://ieeexplore.ieee.org/abstract/document/9304368)
            %%    [31/07/2026] - WGV: the sign of -24.3 for variable x3 in second row of A was corrected.
            %%                   After reviewing, the main paper of 1960 mentioned in You and Dai. The gold
            %%                   standard results are taken also from the paper of 1960 (Land and Doig)
            t = sprintf('%s - nonconvex 7-d MIQCQP : ', names{k});
            if does_nonconv(k)
                H = sparse([1;2],[2;1],[1/2;1/2],7,7,2);
                c = [-77.9 -76.8 -89.6 -97.1 -31.3 0 0]';
                Q = {sparse([1;2;3;3],[3;3;1;2],[-1/2;1/2;-1/2;1/2],7,7,4)};
                B  = zeros(1,7);
                lq = [];
                uq = 0;
                A = [10.9  3.6  -40.8  43.9  7.1   1 0
                     -86.8 32.7 24.3  13.8  -12.6 0 1
                     60.9  68.9  69.0  -56.9 22.5  0 0];
                u = [82.3 77.3 86.5]';
                l = u;
                xmin = zeros(7,1);
                xmax = Inf(7,1);
                x0 = [];
                vtype = 'IIICCCC';
                om = opt_model();
                om.init_set_types();
                om.var.add('x', 7, x0, xmin, xmax, vtype);
                om.qdc.add(om.var, 'cost', H, c);
                om.qcn.add(om.var, 'xQx', Q, B, lq, uq);
                om.lin.add(om.var, 'Ax', A, l, u);
                [x, f, s, out, lam] = om.solve(opt);
                t_is(s, 1, 12, [t 'success']);
                t_is(x([1:6]), [1 0 4 5.0702 1.6930 0]', 4, [t 'x']);
                t_is(f, -981.6023, 4, [t 'f']);

                %% 7) Same previous relaxing integer variable
                t = sprintf('%s - nonconvex 7-d MIQCQP (integer relaxed): ', names{k});
                om.var.set_params('x','v0',ones(7,1));
                [x, f, s, out, lam] = om.solve(opt_r);
                t_is(s, 1, 12, [t 'success']);
                t_is(x, [1.4960 0 5.0210 6.1697 0 0 0]', 4, [t 'x']);
                t_is(f, -1165.50, 2, [t 'f']);

                %% 8) Same problem fixing integer variables to their optimal values
                t = sprintf('%s - nonconvex 7-d MIQCQP (integer fixed) : ', names{k});
                x0 = [1 0 4 zeros(1,4)]';
                om.var.set_params('x','v0',x0);
                [x, f, s, out, lam] = om.solve(opt_f);
                t_is(s, 1, 12, [t 'success']);
                t_is(x([1:6]), [1 0 4 5.0702 1.6930 0]', 4, [t 'x']);
                t_is(f, -981.6023, 4, [t 'f']);
            else
                t_skip(nmiqcqp_nonconvex, sprintf('%s does not handle nonconvex MIQCQP problems', names{k}));
            end
        else
            t_skip(nmiqcqp_convex+nmiqcqp_nonconvex, sprintf('%s does not handle MIQCQP problems', names{k}));
        end
    end
end

t_end;

function TorF = have_miqp_solver()
TorF = have_feature('cplex') || have_feature('gurobi') || have_feature('mosek');
