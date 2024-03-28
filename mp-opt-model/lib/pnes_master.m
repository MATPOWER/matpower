function varargout = pnes_master(fcn, x0, opt)
% pnes_master - Parameterized Nonlinear Equation Solver wrapper function.
% ::
%
%   [X, F, EXITFLAG, OUTPUT, JAC] = PNES_MASTER(FCN, X0, OPT)
%   [X, F, EXITFLAG, OUTPUT, JAC] = PNES_MASTER(PROBLEM)
%   A common wrapper function for numerical continuation methods for
%   solving parameterized nonlinear equations. Traces the solutions of
%   a parameterized nonlinear equation f(x) = 0, beginning from a starting
%   point x0, where f(x) has dimension n and x has dimension n+1.
%
%   In the current implementation, the last element of x is taken to be
%   the parameter lambda, where lambda = 0 corresponds to the base solution.
%
%   Inputs:
%       FCN : handle to function that evaluates the function f(x) to
%           be solved and its Jacobian, J(x). Calling syntax for this
%           function is:
%               f = FCN(x)
%               [f, J] = FCN(x)
%           For a parameterized function, f is n x 1, x is (n+1) x 1,
%           and J is the n x (n+1) matrix of partial derivatives of
%           f (rows) w.r.t. x (cols).
%       X0 : starting value, x0, of vector x ((n+1) x 1)
%       OPT : optional options structure with the following fields,
%           all of which are also optional (default values shown in
%           parentheses)
%           alg ('DEFAULT') : determines which solver to use
%               'DEFAULT' : automatic, currently there is only one
%               solver implementation, a predictor/corrector method
%           verbose (0) - controls level of progress output displayed
%               0 = no progress output
%               1-5 = increasing levels of progress output
%           nleqs_opt - options struct for NLEQS_MASTER used for
%               corrector stage (see NLEQS_MASTER for details), default
%               sets nleqs_opt.verbose to 0, otherwise to 2 if OPT.verbose > 4
%           solve_base (1) : 0/1 flag that determines whether or not to
%               run a corrector stage for initial solution point, x0
%           parameterization (3) - choice of parameterization
%               1 - natural
%               2 - arc len
%               3 - pseudo arc len
%           stop_at ('NOSE') - determines stopping criterion
%               'NOSE'     - stop when limit or nose point is reached
%               'FULL'     - trace full continuation curve
%               <lam_stop> - stop upon reaching specified target lambda value
%           max_it (2000) - maximum number of continuation steps
%           step (0.05) - continuation step size
%           adapt_step (0) - toggle adaptive step size feature
%               0 - adaptive step size disabled
%               1 - adaptive step size enabled
%           adapt_step_damping (0.7) - damping factor for adaptive step sizing
%           adapt_step_tol (1e-3) - tolerance for adaptive step sizing
%           adapt_step_ws (1) - scale factor for default initial step size
%               when warm-starting with adaptive step size enabled
%           step_min (1e-4) - minimum allowed step size
%           step_max (0.2) - maximum allowed step size
%           default_event_tol (1e-3) - default tolerance for event functions
%           target_lam_tol (0) - tolerance for target lambda detection, 0 means
%               use the value of default_event_tol
%           nose_tol (0) - tolerance for nose point detection, 0 means use
%               the value of default_event_tol
%           events (<empty>) - cell array of specs for user-defined event
%               functions, passed as MY_EVENTS arg to PNE_REGISTER_EVENTS
%               (see PNE_REGISTER_EVENTS for details).
%           callbacks (<empty>) - cell array of specs for user-defined callback
%               functions, to be passed as MY_CBACKS arg to
%               PNE_REGISTER_CALLBACKS (see PNE_REGISTER_CALLBACKS for details).
%           output_fcn (<empty>) - handle to custom output function, called by
%               PNE_CALLBACK_DEFAULT
%           plot - options for plotting of continuation curve by
%               PNE_CALLBACK_DEFAULT
%               .level (0) - control plotting of continuation curve
%                   0 - do not plot continuation curve
%                   1 - plot when completed
%                   2 - plot incrementally at each continuation step
%                   3 - same as 2, with 'pause' at each continuation step
%               .idx (<empty>) - index of quantity to plot, passed to yfcn()
%               .idx_default (<empty>) - function to provide default value
%                   for idx, if none provided
%               .xname ('lam') - name of output field holding values that
%                   determine horizontal coordinates of plot
%               .yname ('x') - name of output field holding values that
%                   determine vertical coordinates of plot
%               .xfcn (<empty>) - handle to function that maps a value from
%                   the field of the OUTPUT indicated by value of plot.xname
%                   to a horizontal coordinate for plotting
%               .yfcn (<empty>) - handle to function that maps a value from
%                   the field of the OUTPUT indicated by value of plot.yname
%                   and an index to be applied to that value into a vertical
%                   coordinate for plotting
%               .xlabel ('\lambda') - label for horizontal axis
%               .ylabel ('Variable Value') - label for vertical axis
%               .title ('Value of Variable %d') - plot title used for plot of
%                   single variable, can use %d as placeholder for var index
%               .title2 ('Value of Multiple Variables') - plot title used for
%                   plot of multiple variables
%               .legend ('Variable %d') - legend label, %d can be used as
%                   placeholder for variable index
%           warmstart (<empty>) - struct containing warm-start state, see
%               warmstart field in OUTPUT below for details of expected
%               fields
%       PROBLEM : The inputs can alternatively be supplied in a single
%           PROBLEM struct with fields corresponding to the input arguments
%           described above: fcn, x0, opt
%
%   Outputs (all optional, except X):
%       X : solution vector x
%       F : final function value, f(x)
%       EXITFLAG : exit flag
%           1 = succeeded
%           0 = failed
%       OUTPUT : output struct with the following fields:
%           corrector - output return value from NLEQS_MASTER from final
%               corrector run (see NLEQS_MASTER for details)
%           iterations - N, total number of continuation steps performed
%           events - struct array of size NE of events detected with fields:
%               k - continuation step at which event was located
%               name - name of detected event
%               idx - index(es) of critical elements in corresponding event
%                   function
%               msg - descriptive text detailing the event
%           done_msg - message describing cause of continuation termination
%           steps - (N+1) row vector of stepsizes taken at each continuation
%               step
%           lam_hat - (N+1) row vector of lambda values from prediction steps
%           lam - (N+1) row vector of lambda values from correction steps
%           max_lam - maximum value of parameter lambda (from OUTPUT.lam)
%           warmstart - optional output with information needed for
%               warm-starting an updated continuation problem, with fields:
%               cont_steps - current value of continuation step counter
%               direction - +1 or -1, for tracing of curve in same or
%                   opposite direction, respectively
%               dir_from_jac_eigs - 0/1 flag to indicate whether to use
%                   the sign of the smallest eigenvalue of the Jacobian to
%                   determine the initial direction
%               x - current solution vector
%               z - current tangent vector
%               xp - previous step solution vector
%               zp - previous step tangent vector
%               parm - function handle for current parameterization function
%               default_parm - function handle for default parameterization fcn
%               default_step - default step size
%               events - current event log, same as OUTPUT.events
%               cbs - struct containing user state information for callbacks
%                   see PNES_CALLBACK_DEFAULT for more details
%           (others) - depends on OPT.output_fcn, by default (i.e. with no
%               explicitly provided output function) includes fields:
%                   x_hat - NX x (N+1) matrix of solution values from
%                       prediction steps
%                   x - NX x (N+1) matrix of solution values from correction
%                       steps
%       JAC : final Jacobian matrix, J(x)
%
%   Calling syntax options:
%       [x, f, exitflag, output, jac] = pnes_master(fcn, x0);
%       [x, f, exitflag, output, jac] = pnes_master(fcn, x0, opt);
%       x = pnes_master(problem);
%               where problem is a struct with fields: fcn, x0, opt
%               where opt is optional
%       x = pnes_master(...);
%       [x, f] = pnes_master(...);
%       [x, f, exitflag] = pnes_master(...);
%       [x, f, exitflag, output] = pnes_master(...);
%       [x, f, exitflag, output, jac] = pnes_master(...);
%
%   Example: (based on https://www.chilimath.com/lessons/advanced-algebra/systems-non-linear-equations/)
%       function [f, J] = f1p(x)
%           f = [  x(1)   + x(2) + 6*x(3) - 1;
%                 -x(1)^2 + x(2)          + 5   ];
%           if nargout > 1
%               J = [1 1 6; -2*x(1) 1 0];
%           end
%       end
%       problem = struct( ...
%           'fcn',  @(x)f1p(x), ...
%           'x0',   [-1; 0; 0], ...
%           'opt',  struct('verbose', 2, 'adapt_step', 1, 'step_max', 10) ...
%       );
%       [x, f, exitflag, output, jac] = pnes_master(problem);
%
% See also pne_callback_default, pne_register_callbacks, pne_register_events.

%   MP-Opt-Model
%   Copyright (c) 2013-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell,
%   Shrirang Abhyankar, Argonne National Laboratory,
%   and Alexander Flueck, IIT
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%%----- input argument handling  -----
%% gather inputs
if nargin == 1 && isstruct(fcn) %% problem struct
    p = fcn;
    fcn = p.fcn;
    x0 = p.x0;
    if isfield(p, 'opt'),   opt = p.opt;    else,   opt = struct(); end
else                            %% individual args
    if nargin < 3
        opt = struct();
    end
end

%% default options
if isfield(opt, 'verbose') && opt.verbose > 4
    nleqs_opt_verbose = 2;
else
    nleqs_opt_verbose = 0;
end
dopts = struct( ...
    'alg',              'DEFAULT', ...  %% algorithm
    'verbose',          0, ...
    'nleqs_opt',        struct('verbose', nleqs_opt_verbose), ...
    'solve_base',       1, ...          %% run corrector for initial point
    'parameterization', 3, ...          %% 1 - natural, 2 - arc len, 3 - pseudo arc len
    'stop_at',          'NOSE', ...     %% 'NOSE', 'FULL', <lam_stop>
    'max_it',           2000, ...       %% maximum number of continuation steps
    'step',             0.05, ...       %% continuation step size
    'step_min',         1e-4, ...       %% minimum allowed step size
    'step_max',         0.2, ...        %% maximum allowed step size
    'adapt_step',       0, ...          %% 0/1 toggle, adaptive step size
    'adapt_step_ws',    1, ...          %% warm start inital step size scale factor
    'adapt_step_damping', 0.7, ...      %% adaptive step sizing damping factor
    'adapt_step_tol',   1e-3, ...       %% adaptive step sizing tolerance
    'default_event_tol',1e-3, ...       %% default event function tolerance
    'target_lam_tol',   0, ...          %% event function tolerance for TARGET_LAM event
    'nose_tol',         0, ...          %% event function tolerance for NOSE event
    'events',           {{}}, ...       %% user-defined event detection functions
    'callbacks',        {{}}, ...       %% user-defined callback functions
    'output_fcn',       [], ...         %% custom output fcn, default callback
    'warmstart',        [], ...         %% struct containing warm-start state
    'plot',             struct( ...     %% used by pne_callback_default() for plotting
        'level',        0, ...          %% 0 - no plot, 1 - final, 2 - steps, 3 - steps w/pause
        'idx',          [], ...         %% index of quantity to plot, passed to yfcn()
        'idx_default',  [], ...         %% fcn to provide default value for idx, if none provided
        'xfcn',         [], ...         %% fcn to compute x-coord from data
        'yfcn',         [], ...         %% fcn to compute x-coord from data, idx
        'title',        'Value of Variable %d', ... %% plot title for single var plot
        'title2',       'Value of Multiple Variables', ...  %% plot title for multiple var plot
        'xname',        'lam', ...      %% name of output field holding x vals
        'yname',        'x', ...        %% name of output field holding y vals
        'xlabel',       '\lambda', ...              %% horizontal axis label
        'ylabel',       'Variable Value', ...       %% vertical axis label
        'legend',       'Variable %d' ...           %% legend label
    ) ...
);
opt = nested_struct_copy(dopts, opt);
%% use opt.default_event_tol for NOSE and TARGET_LAM events, unless specified
if opt.target_lam_tol == 0
    opt.target_lam_tol = opt.default_event_tol;
end
if opt.nose_tol == 0
    opt.nose_tol = opt.default_event_tol;
end
if opt.max_it == 0      %% zero means use the default
    opt.max_it == dopts.max_it; 
end

%% initialize
warmstarted = ~isempty(opt.warmstart);
s = struct( ...         %% container struct for various variables, flags
    'done',     0, ...      %% flag indicating continuation has terminated
    'done_msg', '', ...     %% termination message
    'warmstart',[], ...     %% warm start state to return when done (to pass
                            ...%% to subsequent warm-started call to PNES_MASTER)
    'rollback', 0, ...      %% flag to indicate a step must be rolled back
    'events',    [], ...    %% struct array for detected events
    'results',  []  );      %% results struct

%% register event and callback functions
switch opt.stop_at
    case 'NOSE'
        my_events = {{'NOSE', @pne_event_nose, opt.nose_tol}, opt.events{:}};
        my_cbacks = {{@pne_callback_nose, 51}, opt.callbacks{:}};
    case 'FULL'
        my_events = {{'TARGET_LAM', @pne_event_target_lam, opt.target_lam_tol}, opt.events{:}};
        my_cbacks = {{@pne_callback_target_lam, 50}, opt.callbacks{:}};
    otherwise   %% numeric, stop at target lam or nose point, whichever is 1st
        my_events = {{'TARGET_LAM', @pne_event_target_lam, opt.target_lam_tol}, ...
                     {'NOSE', @pne_event_nose, opt.nose_tol}, ...
                        opt.events{:}};
        my_cbacks = {{@pne_callback_target_lam, 50}, ...
                     {@pne_callback_nose, 51}, ...
                        opt.callbacks{:}};
end
my_cbacks{end+1} = {@pne_callback_default, 0};
reg_ev = pne_register_events(my_events, opt);   %% registered event functions
reg_cb = pne_register_callbacks(my_cbacks);     %% registered callback functions
nef = length(reg_ev);   %% number of registered event functions
ncb = length(reg_cb);   %% number of registered callback functions

%% initialize continuation step counter
if warmstarted
    cont_steps = opt.warmstart.cont_steps + 1;
    if opt.verbose
        fprintf('... CONTINUATION RESUMED\n');
    end
else
    cont_steps = 0;
    if opt.verbose
        v = mpomver('all');
        fprintf('\nMP-Opt-Model Version %s, %s', v.Version, v.Date);
        fprintf(' -- Predictor/Corrector Continuation Method\n');
    end
end

%% solve corrector step for base point
if opt.solve_base && ~warmstarted
    cfcn = @(xx)pne_corrector_fcn(xx, fcn, @pne_pfcn_natural, x0, 0, []);
    [x, f, exitflag, out] = nleqs_master(cfcn, x0, opt.nleqs_opt);
    if exitflag
        if opt.verbose > 1
            fprintf('step %3d  :                          lambda = %6.3f, %2d corrector steps\n', cont_steps, x0(end), out.iterations);
        end
    else
        s.done = 1;
        s.done_msg = sprintf('base solution did not converge in %d iterations', out.iterations);
        if opt.verbose
            fprintf('%s\n', s.done_msg);
        end
    end
else
    x = x0;     %% ignored for warmstart, overwritten by warmstart.cx
end

%% initialize numerical continuation
if ~s.done
    locating = 0;   %% flag to indicate that an event interval was detected,
                    %% but the event has not yet been located
    rb_cnt_ef = 0;  %% counter for rollback steps triggered by event function intervals
    rb_cnt_cb = 0;  %% counter for rollback steps triggered directly by callbacks

    if warmstarted
        manual_direction_switch = 0;    %% set to 1 to manually prompt for
                                        %% direction change upon warmstart
        ws = opt.warmstart;
        x = ws.x;           %% starting value solution vector
        z = ws.z;           %% starting value of tangent vector
        step = 0;           %% re-solve current point
        parm = ws.parm;
        direction = ws.direction;
        default_parm = ws.default_parm;
        default_step = ws.default_step;
        cbs = ws.cbs;
        event_log = ws.events;
        if manual_direction_switch
            %% decide whether to switch directions
            reply = input('Switch directions? Y/N [N]:','s');
            if strcmp(upper(reply), 'Y')
                direction = -direction;
            end
        elseif isfield(ws, 'dir_from_jac_eigs') && ws.dir_from_jac_eigs
            %% attempt to determine direction from smalles Jacobian eigenvalue
            [~, J] = fcn(x);
            eigs_opt.tol = 1e-3;
            eigs_opt.maxit = 2*length(x);
            direction = sign(z(end) * ...
                        min(real(eigs(J(:,1:end-1), 1, 'sr', eigs_opt))));
        end

        if opt.adapt_step   %% hey, maybe slow down, things might have changed
            default_step = default_step * opt.adapt_step_ws;
        end
    else
        %% initialize parameterization function
        switch opt.parameterization
            case 1
                parm = @pne_pfcn_natural;           %% NAT
            case 2
                parm = @pne_pfcn_arc_len;           %% ARC
            case 3
                parm = @pne_pfcn_pseudo_arc_len;    %% PAL
            otherwise
                error('pnes_master: OPT.parameterization (= %d) must be 1, 2, or 3', opt.parameterization);
        end

        %% finish initializing tangent vector
        direction = 1;
        z0 = zeros(length(x), 1); z0(end) = direction;  %% +ve lambda direction
        z = pne_tangent(x, x, z0, fcn, parm, direction);

        step = opt.step;
        default_step = step;
        default_parm = parm;
        cbs = [];
        event_log = [];
    end

    %% initialize state struct for current continuation step
    cx = struct( ...        %% current state
        'x_hat',        x, ...      %% predicted solution value
        'x',            x, ...      %% corrected solution value
        'z',            z, ...      %% normalized tangent vector
        'default_step', default_step, ...   %% default step size
        'default_parm', default_parm, ...   %% default parameterization
        'this_step', [], ...        %% step size for this step only
        'this_parm', [], ...        %% parameterization for this step only
        'step', step, ...           %% current step size
        'parm', parm, ...           %% current parameterization
        'events', event_log, ...    %% event log
        'cbs', cbs, ...             %% user-defined callback state
        'efv', [] ...               %% event function values
    );

    %% initialize event function values
    cx.efv = cell(nef, 1);
    for k = 1:nef
        cx.efv{k} = reg_ev(k).fcn(cx, opt);
    end

    if warmstarted  %% no need to initialize callbacks
        %% initialize state for previous continuation step
        px = cx;
        px.x = ws.xp;   %% use warm start value for solution value
        px.z = ws.zp;   %% use warm start value for tangent
    else
        %% invoke callbacks - "initialize" context
        for k = 1:ncb
            [nx, cx, s] = reg_cb(k).fcn(cont_steps, cx, cx, cx, s, opt);
        end
        cont_steps = cont_steps + 1;

        %% check for case with base and target the same
        %% evaluate function at lambda = 0 (base)
        if opt.solve_base
            fb = f(1:end-1);
        else
            fb = fcn(x);
            exitflag = 1;
        end

        %% evaluate function at lambda = 1 (target)
        xt = x;
        xt(end) = 1;
        ft = fcn(xt);
        if norm(fb - ft, Inf) < 1e-12
            s.done = 1;
            s.done_msg = 'base and target functions are identical';
        end

        %% initialize state for previous continuation step
        px = cx;
    end
end

%%-----  run numerical continuation  -----
while ~s.done
    %% initialize next candidate with current state
    nx = cx;

    %% predictor step
    nx.x_hat = cx.x + cx.step * cx.z;

    %% corrector step
    cfcn = @(xx)pne_corrector_fcn(xx, fcn, cx.parm, cx.x, cx.step, cx.z);
    [nx.x, f, exitflag, out] = nleqs_master(cfcn, nx.x_hat, opt.nleqs_opt);
    if ~exitflag        %% corrector failed
        s.done = 1;
        s.done_msg = sprintf('Corrector did not converge in %d iterations.', out.iterations);
        if opt.verbose
            fprintf('step %3d  : %s stepsize = %-9.3g lambda = %6.3f  corrector did not converge in %d iterations\n', cont_steps, pne_ptag(cx.parm), cx.step, nx.x(end), out.iterations);
        end
        cont_steps = max(cont_steps - 1, 1);    %% go back to last step, but not to 0
        break;
    end

    %% compute new tangent direction, based on current or prev state: tx
    if nx.step == 0     %% if this is a re-do step, cx and nx are the same
        tx = px;            %% so use px as the previous state
    else                %% otherwise
        tx = cx;            %% use cx as the previous state
    end
    nx.z = pne_tangent(nx.x, tx.x, tx.z, fcn, nx.parm, direction);
    direction = 1;      %% continue in same direction

    %% detect events
    for k = 1:nef
        nx.efv{k} = reg_ev(k).fcn(nx, opt); %% update event function values
    end
    [s.rollback, s.events, nx.efv] = ...
        pne_detect_events(reg_ev, nx.efv, cx.efv, nx.step);

    %% adjust step-size to locate event function zero, if necessary
    if s.rollback               %% current step overshot
        %% roll back & initialize next step size based on rollback and previous
        rx = nx;                    %% save state we're rolling back from
        rx_evnts = s.events;        %% and critical event info
        cx.this_step = s.events.step_scale * rx.step;
        cx.this_parm = rx.parm;     %% keep same parameterization as last step
        locating = 1;               %% enter "locating" mode (or stay in it)
        rb_cnt_ef = rb_cnt_ef + 1;  %% increment rollback counter for ef intervals
        if rb_cnt_ef > 26
            s.done = 1;
            s.done_msg = sprintf('Could not locate %s event!', s.events.name);
        end
        if opt.verbose > 3
            loc_msg = sprintf('OVERSHOOT  : f = [%g, <<%g>>], step <-- %.4g', ...
                        cx.efv{s.events.eidx}(s.events.idx(1)), ...
                        rx.efv{s.events.eidx}(s.events.idx(1)), cx.this_step);
        end
    elseif locating
        if s.events(1).zero      %% found the zero!
            %% step size will be reset to previously used default step size
            locating = 0;           %% exit "locating" mode
            rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
            if opt.verbose > 3
                loc_msg = sprintf('ZERO!      : f = %g, step <-- %.4g', ...
                    nx.efv{rx_evnts.eidx}(rx_evnts.idx(1)), nx.default_step);
            end
        else                    %% prev rollback undershot
            %% initialize next step size based on critical event function
            %% values from prev rollback step and current step
            rx_efv = rx.efv{rx_evnts.eidx}(rx_evnts.idx(1));
            cx_efv = nx.efv{rx_evnts.eidx}(rx_evnts.idx(1));
            step_scale = cx_efv / (cx_efv - rx_efv);
            nx.this_step = step_scale * (rx.step - nx.step);
            rb_cnt_ef = 0;          %% reset rollback counter for ef intervals
            if opt.verbose > 3
                loc_msg = sprintf('UNDERSHOOT : f [<<%g>>, %g], step <-- %.4g', ...
                    cx_efv, rx_efv, nx.this_step);
            end
        end
    else                    %% normal step, not locating anything
        loc_msg = '';
    end

    %% invoke callbacks - "iterations" context
    rb = s.rollback;
    for k = 1:ncb
        [nx, cx, s] = reg_cb(k).fcn(cont_steps, nx, cx, px, s, opt);
    end
    if ~rb && s.rollback    %% rollback triggered by callback (vs event function interval)
        rb_cnt_cb = rb_cnt_cb + 1;  %% increment rollback counter for callbacks
        if rb_cnt_cb > 26
            s.done = 1;
            s.done_msg = 'Too many rollback steps triggered by callbacks!';
        end
    else
        rb_cnt_cb = 0;              %% reset rollback counter for callbacks
    end

    %% print iteration information
    if opt.verbose > 1
        %% set label for rollback step counter
        if rb_cnt_ef
            sub_step = char('a' + rb_cnt_ef - 1);
        elseif rb_cnt_cb
            sub_step = char('A' + rb_cnt_cb - 1);
        else
            sub_step = ' ';
        end

        fprintf('step %3d%s : %s stepsize = %-9.3g lambda = %6.3f', cont_steps, sub_step, pne_ptag(cx.parm), cx.step, nx.x(end));
        if opt.verbose < 5
            fprintf('  %2d corrector steps', out.iterations);
        end
        if s.rollback
            fprintf(' ^ ROLLBACK\n');
        else
            fprintf('\n');
        end
        if opt.verbose > 3 && ~isempty(loc_msg)
            fprintf('    LOCATING -- %s\n', loc_msg);
        end
    end

    %% log events
    for k = 1:length(s.events)
        if s.events(k).log
            e = struct( 'k', cont_steps, ...
                        'name', s.events(k).name, ...
                        'idx', s.events(k).idx, ...
                        'msg', s.events(k).msg   );
            if isempty(nx.events)
                nx.events = e;
            else
                nx.events(end+1) = e;
            end
        end
        if (opt.verbose > 2 && s.events(k).log) || ...
                (opt.verbose > 3 && s.events(k).eidx)
            fprintf('    %s\n', s.events(k).msg);
        end
    end

    %% adapt stepsize if requested and not terminating, locating a zero
    %% or warm starting
    if opt.adapt_step && ~s.done && ~locating && ~s.events(1).zero && nx.step ~= 0
        pred_error = norm(nx.x - nx.x_hat, Inf);

        %% new nominal step size is current size * tol/err, but we reduce
        %% the change from the current size by a damping factor and limit
        %% increases to a factor of 2
        step_scale = min(2, 1 + opt.adapt_step_damping * ...
                        (opt.adapt_step_tol/pred_error - 1));
        nx.default_step = nx.step * step_scale;

        %% limit step-size
        if nx.default_step > opt.step_max
            nx.default_step = opt.step_max;
        end
        if nx.default_step < opt.step_min
            nx.default_step = opt.step_min;
        end
    end

    %% if this is a normal step
    if ~s.rollback
        px = cx;    %% save current state before update
        cx = nx;    %% update current state to next candidate
        if ~s.done
            if cont_steps >= opt.max_it
                s.done = 1;
                s.done_msg = sprintf('Reached maximun number of continuation steps (opt.max_it = %d)', opt.max_it);
            else
                cont_steps = cont_steps + 1;
            end
        end
    end

    %% set step size and parameterization, from one-time or defaults
    if isempty(cx.this_step)
        cx.step = cx.default_step;
    else
        cx.step = cx.this_step;
        cx.this_step = [];      %% disable for next time
    end
    if isempty(cx.this_parm)
        cx.parm = cx.default_parm;
    else
        cx.parm = cx.this_parm;
        cx.this_parm = [];      %% disable for next time
    end
end     %% while ~s.done

%% invoke callbacks - "final" context
s.results = struct();   %% initialize results struct
for k = 1:ncb
    [nx, cx, s] = reg_cb(k).fcn(-cont_steps, nx, cx, px, s, opt);
end
output = s.results;
output.done_msg = s.done_msg;
output.events = cx.events;  %% copy eventlog to results
output.corrector = out;     %% output from last corrector run

%% prepare to exit
if isempty(s.warmstart)
    if opt.verbose
        fprintf('CONTINUATION TERMINATION: %s\n', s.done_msg);
    end
else
    %% save warmstart values
    ws = s.warmstart;
    ws.cont_steps = cont_steps;
    ws.direction = direction;

    %% from state at current step
    ws.x = cx.x;            %% state from current step
    ws.z = cx.z;            %% tangent vector from current step
    ws.parm = cx.parm;
    ws.default_parm = cx.default_parm;
    ws.default_step = cx.default_step;
    ws.events = cx.events;
    ws.cbs = cx.cbs;

    %% from state at previous step
    ws.xp = px.x;           %% state from previous step
    ws.zp = px.z;           %% tangent vector from previous step

    output.warmstart = ws;

    if opt.verbose
        fprintf('%s : CONTINUATION SUSPENDED ...\n', s.done_msg);
    end
end

%% output arguments
if nargout > 4
    [f, J] = fcn(cx.x);
elseif nargout > 1
    f = fcn(cx.x);
end
varargout{1} = cx.x;
if nargout > 1
    varargout{2} = f;
    if nargout > 2
        varargout{3} = exitflag;
        if nargout > 3
            varargout{4} = output;
            if nargout > 4
                varargout{5} = J;
            end
        end
    end
end


%%-----  pne_corrector_fcn  -----
%% fcn(x) combined with parameterization constraint
function [fp, dfp] = pne_corrector_fcn(x, fcn, parm, cx_x, step, z)
if nargout < 2
    fp = [ fcn(x); parm(x, cx_x, step, z) ];
else
    [f, df] = fcn(x);
    [p, dp] = parm(x, cx_x, step, z);
    fp = [f; p];
    dfp = [df; dp];
end


%%-----  pne_tangent  -----
%% find normalized tangent vector
function z = pne_tangent(x, xp, zp, fcn, parm, direction)
[f, df] = fcn(x);
[p, dp] = parm(x, xp, 0, zp);
rhs = [ zeros(length(f), 1); direction ];
z = [df; dp] \ rhs;
z = z / norm(z);    %% normalize it


%%-----  pne_ptag  -----
%% return 3-letter string to indicate parameterization scheme
%% NAT - natural, ARC - arc length, PAL - pseudo arc length
function ptag = pne_ptag(parm)
switch func2str(parm)
    case 'pne_pfcn_natural'
        ptag = 'NAT';
    case 'pne_pfcn_arc_len'
        ptag = 'ARC';
    case 'pne_pfcn_pseudo_arc_len'
        ptag = 'PAL';
end
