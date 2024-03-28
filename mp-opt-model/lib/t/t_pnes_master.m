function t_pnes_master(quiet)
% t_pnes_master - Tests of PNE solvers via pnes_master.

%   MP-Opt-Model
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

if nargin < 1
    quiet = 0;
end

%%  alg         name        check       opts
cfg = {
    {'DEFAULT', 'default',  []          []  },
};

n = 188;

t_begin(n*length(cfg), quiet);

for k = 1:length(cfg)
    alg   = cfg{k}{1};
    name  = cfg{k}{2};
    check = cfg{k}{3};
    opts  = cfg{k}{4};
    if ~isempty(check) && ~have_feature(check)
        t_skip(n, sprintf('%s not installed', name));
    else
        %% default start point
        x0 = [-1;0;0];
        opt = struct( ...
            'verbose', 0, ...
            'alg', alg, ...
            'stop_at', 0.7, ...
            'parameterization', 1, ...
            'nleqs_opt', struct('tol', 1e-9), ...
            'adapt_step', 0, ...
            'step_max', 10, ...
            'target_lam_tol', 1e-6, ...
            'nose_tol', 1e-6, ...
            'adapt_step_tol', 1e-2 );
%         opt.plot = struct('level', 2, 'idx', 2);

        t = sprintf('%s - TARGET_LAM = 0.7 (natural) : ', name);
        [x, f, e, out, jac] = pnes_master(@f1p, x0, opt);
        it = 14;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 0.7, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'TARGET_LAM', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it), [t 'out.done_msg']);
        t_ok(isfield(out, 'corrector'), [t 'out.corrector exists']);
        t_is(out.corrector.iterations, 3, 12, [t 'out.corrector.iterations']);
        t_str_match(out.corrector.message, sprintf('Newton''s method converged in %d iterations.', 3), [t 'out.corrector.message']);

        t = sprintf('%s - TARGET_LAM = 0.7 (arc len) : ', name);
        opt.adapt_step = 1;
        opt.parameterization = 2;
        [x, f, e, out, jac] = pnes_master(@f1p, x0, opt);
        it = 9;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 0.7, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'TARGET_LAM', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it), [t 'out.done_msg']);

        t = sprintf('%s - TARGET_LAM = 0.7 (pseudo arc len) : ', name);
        opt.parameterization = 3;
        [x, f, e, out, jac] = pnes_master(@f1p, x0, opt);
        it = 9;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-1.931782106; -1.268217894; 0.7], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 0.7, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'TARGET_LAM', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached desired lambda 0.7 in %d continuation steps', it), [t 'out.done_msg']);

        t = sprintf('%s - FULL : ', name);
        opt.stop_at = 'FULL';
        p = struct('fcn', @f1p, 'x0', x0, 'opt', opt);
        [x, f, e, out, jac] = pnes_master(p);
        it = 34;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [2;-1;0], 10, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 1.04127275, 8, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'TARGET_LAM', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Traced full continuation curve in %d continuation steps', it), [t 'out.done_msg']);

        t = sprintf('%s - FULL (max_it) : ', name);
        opt.max_it = 25;
        p = struct('fcn', @f1p, 'x0', x0, 'opt', opt);
        [x, f, e, out, jac] = pnes_master(p);
        opt.max_it = 2000;
        it = 25;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [0.64152370; -4.588447338; 0.82448727], 8, [t 'x - final']);
        t_is(f, [0;0], 8, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 1.04127275, 8, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 0, 12, [t 'length(out.events)']);
        t_str_match(out.done_msg, sprintf('Reached maximun number of continuation steps (opt.max_it = %d)', it), [t 'out.done_msg']);

        t = sprintf('%s - NOSE (arc len): ', name);
        p.opt.stop_at = 'NOSE';
        p.opt.parameterization = 2;
        [x, f, e, out, jac] = pnes_master(p);
        it = 18;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 1.04166666667, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'NOSE', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it), [t 'out.done_msg']);

        t = sprintf('%s - NOSE (pseudo arc len) : ', name);
        p.opt.parameterization = 3;
        [x, f, e, out, jac] = pnes_master(p);
        it = 18;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [-3;4;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 1.04166666667, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'NOSE', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it), [t 'out.done_msg']);

        t = sprintf('%s - NOSE (opp dir) : ', name);
        p.x0 = [1;-1;0];
        [x, f, e, out, jac] = pnes_master(p);
        it = 20;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-0.5; -4.75; 1.04166667], 5, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_is(out.x(:,1), [2;-1;0], 8, [t 'out.x(:,1)']);
        t_is(out.max_lam, 1.04166666667, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.x), [3,it+1], 12, [t 'size(out.x)']);
        t_is(size(out.x_hat), [3,it+1], 12, [t 'size(out.x_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 1, 12, [t 'length(out.events)']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'NOSE', [t 'out.events.name']);
        t_str_match(out.done_msg, sprintf('Reached limit in %d continuation steps, lambda = 1.042.', it), [t 'out.done_msg']);

        t = sprintf('%s - FULL warmstart, before SWITCH : ', name);
        opt.parameterization = 3;
        opt.adapt_step = 1;
        opt.adapt_step_ws = 0.25;
        opt.stop_at = 'FULL';
        opt.adapt_step_tol = 5e-3;
        opt.callbacks = {@pne_callback_test1};
        opt.events = {{'SWITCH!', @pne_event_test1, 1e-6}};
        opt.output_fcn = @pne_output_fcn_test1;
        % opt.verbose = 4;
        % opt.plot = struct('level', 2, 'yname', 'y', 'idx', [1;2]);
        x0 = [-1;0;0];
        [x, f, e, out, jac] = pnes_master(@f1p, x0, opt);
        it = 10;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [-2; -1; 2/3], 6, [t 'x - final']);
        t_is(f, [0;0], 10, [t 'f']);
        t_ok(isfield(out, 'warmstart') && ~isempty(out.warmstart), [t 'out.warmstart exists']);
        t_is(out.max_lam, 2/3, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.y), [2,it+1], 12, [t 'size(out.y)']);
        t_is(size(out.y_hat), [2,it+1], 12, [t 'size(out.y_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_ok(isstruct(out.events), [t 'out.events is struct']);
        t_is(out.events.k, it, 12, [t 'out.events.k']);
        t_is(out.events.idx, 1, 12, [t 'out.events.idx']);
        t_str_match(out.events.name, 'SWITCH!', [t 'out.events.name']);
        t_str_match(out.events.msg, 'ZERO detected for SWITCH! event', [t 'out.events.msg']);
        t_is(out.warmstart.cont_steps, it, 12, [t 'out.warmstart.cont_steps']);
        t_is(out.warmstart.default_step, 0.64432407, 8, [t 'out.warmstart.default_step']);
        t_ok(isa(out.warmstart.parm, 'function_handle'), [t 'out.warmstart.parm is function']);
        t_ok(isa(out.warmstart.default_parm, 'function_handle'), [t 'out.warmstart.default_parm is function']);
        t_ok(isstruct(out.warmstart.cbs), [t 'out.warmstart.cbs is struct']);
        t_ok(isstruct(out.warmstart.cbs.default), [t 'out.warmstart.cbs.default is struct']);
        t_is(out.warmstart.cbs.default.iterations, it, 12, [t 'out.warmstart.cbs.default.iterations']);
        t_ok(isstruct(out.warmstart.events), [t 'out.warmstart.events is struct']);
        t_is(out.warmstart.events.k, it, 12, [t 'out.warmstart.events.k']);
        t_is(out.warmstart.events.idx, 1, 12, [t 'out.warmstart.events.idx']);
        t_str_match(out.warmstart.events.name, 'SWITCH!', [t 'out.warmstart.events.name']);
        t_str_match(out.warmstart.events.msg, 'ZERO detected for SWITCH! event', [t 'out.warmstart.events.msg']);
        t_str_match(out.done_msg, sprintf('Reached switching point in %d continuation steps', it), [t 'out.done_msg']);
        t_ok(isfield(out, 'corrector'), [t 'out.corrector exists']);
        t_is(out.corrector.iterations, 3, 12, [t 'out.corrector.iterations']);
        t_str_match(out.corrector.message, sprintf('Newton''s method converged in %d iterations.', 3), [t 'out.corrector.message']);

        t = sprintf('%s - FULL warmstart, after SWITCH : ', name);
        ws = out.warmstart;
        ws.x  = [ws.x(end)/2;  ws.x];
        ws.xp = [ws.xp(end)/2; ws.xp];
        ws.z = [0; ws.z];
        ws.zp = [0; ws.zp];
        x0 = ws.x;      %% ignored for warm start
        opt.warmstart = ws;
        opt.events = {{'SNOUT!', @pne_event_test2, 1e-6}};
        opt.callbacks = {};
        [x, f, e, out, jac] = pnes_master(@f2p, x0, opt);
        it = 49;
        t_is(e, 1, 12, [t 'exitflag']);
        t_is(x, [0; 6; -5; 0], 6, [t 'x - final']);
        t_is(f, [0;0;0], 10, [t 'f']);
        t_is(out.y(:,1), [-3;4], 8, [t 'out.y(:,1)']);
        t_is(out.max_lam, 1.04166666667, 10, [t 'out.max_lam']);
        t_is(out.iterations, it, 12, [t 'out.iterations']);
        t_is(size(out.lam), [1,it+1], 12, [t 'size(out.lam)']);
        t_is(size(out.lam_hat), [1,it+1], 12, [t 'size(out.lam_hat)']);
        t_is(size(out.y), [2,it+1], 12, [t 'size(out.y)']);
        t_is(size(out.y_hat), [2,it+1], 12, [t 'size(out.y_hat)']);
        t_is(size(out.steps), [1,it+1], 12, [t 'size(out.steps)']);
        t_is(length(out.events), 3, 12, [t 'length(out.events)']);
        t_is(out.events(1).k, 10, 12, [t 'out.events(1).k']);
        t_is(out.events(1).idx, 1, 12, [t 'out.events(1).idx']);
        t_str_match(out.events(1).name, 'SWITCH!', [t 'out.events(1).name']);
        t_str_match(out.events(1).msg, 'ZERO detected for SWITCH! event', [t 'out.events(1).msg']);
        t_is(out.events(2).k, 28, 12, [t 'out.events(2).k']);
        t_is(out.events(2).idx, 1, 12, [t 'out.events(2).idx']);
        t_str_match(out.events(2).name, 'SNOUT!', [t 'out.events(2).name']);
        t_str_match(out.events(2).msg, 'ZERO detected for SNOUT! event', [t 'out.events(2).msg']);
        t_is(out.events(3).k, it, 12, [t 'out.events(3).k']);
        t_is(out.events(3).idx, 1, 12, [t 'out.events(3).idx']);
        t_str_match(out.events(3).name, 'TARGET_LAM', [t 'out.events(3).name']);
        t_str_match(out.events(3).msg, 'ZERO detected for TARGET_LAM event', [t 'out.events(3).msg']);
        t_str_match(out.done_msg, sprintf('Traced full continuation curve in %d continuation steps', it), [t 'out.done_msg']);
        t_ok(isfield(out, 'corrector'), [t 'out.corrector exists']);
        t_is(out.corrector.iterations, 3, 12, [t 'out.corrector.iterations']);
        t_str_match(out.corrector.message, sprintf('Newton''s method converged in %d iterations.', 3), [t 'out.corrector.message']);
    end
end

t_end;

% lam = 1;
% m = 6;
% xs = [-m:0.2:m];
% ys = xs;
% [xx, yy] = meshgrid(xs, ys);
% zz0 = zeros([size(xx), 2]);
% zz = zz0;
% for i = 1:length(xs)
%     for j = 1:length(ys)
%         zz(i, j, :) = f1p([xx(i, j); yy(i, j); lam]);
% %         zz(i, j, :) = f2([xx(i, j); yy(i, j)]);
%     end
% end
% 
% figure
% ax = gca;
% surf(ax, xx, yy, squeeze(zz(:, :, 1)))
% hold on
% surf(ax, xx, yy, squeeze(zz(:, :, 2)))
% surf(ax, m*[-1 -1; 1 1], m*[-1 1; -1 1], [0 0; 0 0])
% hold off


%% 2-d problem with 2 solutions
%% from https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1(x)
f = [  x(1)   + x(2) - 1;
      -x(1)^2 + x(2) + 5    ];
if nargout > 1
    J = [1 1; -2*x(1) 1];
end

%% parameterized 2-d problem
%% based on https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/
function [f, J] = f1p(x)
if nargout < 2
    f = f1(x(1:2));
else
    [f, J] = f1(x(1:2));
    J = [J [6;0]];
end
f = f + [6*x(3); 0];

%% parameterized 3-d problem, based on f1p
function [f, J] = f2p(x)
x1 = [x(3)+2; x(2)-2; x(4)];
if nargout < 2
    f = [2*x(1)-x(4); f1p(x1)];
else
    [f, JJ] = f1p(x1);
    J = [  2      0        0       -1;
          [0; 0] JJ(:, 2) JJ(:, 1) JJ(:, 3) ];
end
f = [ 2*x(1)-x(4); f ];

%% example custom event function 1 (target lambda == 2/3)
function efv = pne_event_test1(cx, opt)
tlam = 2/3;
efv = cx.x(end) - tlam;

%% example custom event function 2 (nose point)
function efv = pne_event_test2(cx, opt)
efv = cx.z(end);

%% example custom callback function (exit for warmstart on SWITCH! event)
function [nx, cx, s] = pne_callback_test1(k, nx, cx, px, s, opt)
if k <= 0 || s.done, return; end    %% skip if initialize, finalize or done
tlam = 2/3;
ev = pne_detected_event(s.events, 'SWITCH!');
if ~isempty(ev)
    if ev.zero              %% prepare to terminate
        s.done = 1;
        s.done_msg = sprintf('Reached switching point in %d continuation steps', k);
        s.warmstart = struct(); %% signal that we want to exit, then resume
    else                    %% set step-size & parameterization to terminate next time
        cx.this_parm = @pne_pfcn_natural;   %% change to natural parameterization
        cx.this_step = tlam - cx.x(end);
        ev.msg = sprintf('%s\n  step %d to overshoot lambda = %g, reduce step size and set natural param', ev.msg, k, tlam);
    end
end

function [names, vals] = pne_output_fcn_test1(x, x_hat)
%% [names, vals] = pne_output_fcn_default(x, x_hat)
%% names = pne_output_fcn_default()
names = {'y_hat', 'y'};
if nargin
    if length(x) == 3
        k = [1;2];
    else %% == 3
        k = [2;3];
    end
    vals = {x_hat(k), x(k)};
end


% %% another 2-d problem
% %% from Christi Patton Luks, https://www.youtube.com/watch?v=pJG4yhtgerg
% function [f, J] = f2(x)
% m = 10;
% m1 = 0; %2;
% m2 = 0; %5;
% f = [  x(1)^2 +   x(1)*x(2)   - 10 + m1;
%        (x(2)   + 3*x(1)*x(2)^2 - 57) / m  ] + m2;
% if nargout > 1
%     J = sparse([    2*x(1)+x(2)    x(1);
%                     (3*x(2)^2) / m (6*x(1)*x(2)+1) / m    ]);
% end
