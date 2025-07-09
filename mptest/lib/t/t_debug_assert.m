function t_debug_assert(quiet)
% t_debug_assert - Test debug_assert function.
% ::
%
%   t_debug_assert
%   t_debug_assert(quiet)

%   Copyright (c) 2024, Ray Zimmerman

if nargin < 1
    quiet = 0;
end

n_tests = 8;

t_begin(n_tests, quiet);

%% save any pre-existing state of DEBUG_MODE global
global DEBUG_MODE;
if ~isempty('DEBUG_MODE')
    need_to_restore = true;
    saved_debug_mode = DEBUG_MODE;
    DEBUG_MODE = [];
    clear('DEBUG_MODE');
else
    need_to_restore = false;
end

t = 'DEBUG_MODE undefined : ';
try
    debug_assert(true, 'one');
    t_ok(true, [t 'true']);
catch
    t_ok(false, [t 'true']);
end
try
    debug_assert(false, 'one');
    t_ok(true, [t 'false']);
catch
    t_ok(false, [t 'false']);
end

t = 'DEBUG_MODE (empty) : ';
global DEBUG_MODE;
try
    debug_assert(true, 'one');
    t_ok(true, [t 'true']);
catch
    t_ok(false, [t 'true']);
end
try
    debug_assert(false, 'one');
    t_ok(true, [t 'false']);
catch
    t_ok(false, [t 'false']);
end

t = 'DEBUG_MODE = false : ';
if quiet
    toggle_debug_mode(false, quiet);
else
    toggle_debug_mode(false);
end
try
    debug_assert(true, 'one');
    t_ok(true, [t 'true']);
catch
    t_ok(false, [t 'true']);
end
try
    debug_assert(false, 'one');
    t_ok(true, [t 'false']);
catch
    t_ok(false, [t 'false']);
end

t = 'DEBUG_MODE = true : ';
if quiet
    toggle_debug_mode([], quiet);
else
    toggle_debug_mode();
end
try
    debug_assert(true, 'one');
    t_ok(true, [t 'true']);
catch
    t_ok(false, [t 'true']);
end
try
    debug_assert(false, 'one');
    t_ok(false, [t 'false']);
catch
    t_ok(true, [t 'false']);
end

t_end;

%% restore any pre-existing state of DEBUG_MODE global
if need_to_restore
    DEBUG_MODE = saved_debug_mode;
else
    DEBUG_MODE = [];
    clear('DEBUG_MODE');
end
