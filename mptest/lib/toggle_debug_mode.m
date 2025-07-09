function toggle_debug_mode(on_off, quiet)
% toggle_debug_mode - Toggles or sets the value of ``DEBUG_MODE`` global var.
% ::
%
%   toggle_debug_mode
%   toggle_debug_mode(true)
%   toggle_debug_mode(false)
%   toggle_debug_mode([], quiet)
%   toggle_debug_mode(true, quiet)
%   toggle_debug_mode(false, quiet)
%
% Used to control whether or not debug_assert will call :func:`assert`
% or do nothing. Displays the message ``DEBUG_MODE is ON`` or
% ``DEBUG_MODE is ON``, depending on the new state, unless the ``quiet``
% argument is present and true.
%
% See also debug_assert.

%   Copyright (c) 2024-2025, Ray Zimmerman

global DEBUG_MODE
if nargin && ~isempty(on_off)
    if on_off
        DEBUG_MODE = true;
    else
        DEBUG_MODE = false;
    end
elseif DEBUG_MODE
    DEBUG_MODE = false;
else
    DEBUG_MODE = true;
end
if nargin < 2 || ~quiet
    if DEBUG_MODE
        fprintf('DEBUG_MODE is ON\n');
    else
        fprintf('DEBUG_MODE is OFF\n');
    end
end
