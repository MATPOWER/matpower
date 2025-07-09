function debug_assert(cond, varargin)
% debug_assert - Calls  :func:`assert` if and only if ``DEBUG_MODE`` is true.
% ::
%
%   debug_assert(cond)
%   debug_assert(cond, msg, A)
%   debug_assert(cond, errID, msg)
%   debug_assert(cond, errID, msg, A)
%
% Calls :func:`assert` if global variable ``DEBUG_MODE`` exists and is
% true, otherwise skips the call to :func:`assert` and does nothing.
%
% See also toggle_debug_mode, assert.

%   Copyright (c) 2024-2025, Ray Zimmerman

global DEBUG_MODE
if DEBUG_MODE
    assert(cond, varargin{:});
end
