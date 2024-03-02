function rv = have_feature(tag, rtype)
% have_feature  Test for optional functionality / version info.
% ::
%
%   TorF = have_feature(tag)
%   TorF = have_feature(tag, toggle)
%   ver_str = have_feature(tag, 'vstr')
%   ver_num = have_feature(tag, 'vnum')
%   rdate = have_feature(tag, 'date')
%   all_info = have_feature(tag, 'all')
%   have_feature(tag, 'clear_cache')
%   have_feature('all', 'clear_cache')
%
% Returns availability, version and release information for optional
% functionality. All information is cached, and the cached values
% returned on subsequent calls. If the functionality exists, an attempt is
% made to determine the release date and version number.
%
% Inputs:
%   tag (char array) : name to identify optional functionality
%   toggle (-1, 0, 1) : numeric second argument, to turn optional functionality on or off
%
%       ==========  ============  ========================================
%       ``toggle``  return type   description and return value
%       ==========  ============  ========================================
%       *<none>*    *boolean*     check availability, return true if available, false if not
%       0           *boolean*     **turn off** optional functionality, return false
%       1           *boolean*     **turn on** optional functionality, return
%                                 true if available, false if not
%       -1          *boolean*     **toggle** on/off state of optional
%                                 functionality, return true if available,
%                                 false if not
%       ==========  ============  ========================================
%   rtype (char array) : string valued second argument determines which
%       value is returned for available functionality, as follows:
%
%       =================  ============  ========================================
%       ``rtype``          return type   description and return value
%       =================  ============  ========================================
%       ``'vstr'``         *char array*  version number as a string
%       ``'vnum'``         *double*      version number as numeric value
%       ``'date'``         *char array*  release date
%       ``'all'``          *struct*      struct with all info (availablity, version numbers, release date)
%       ``'clear_cache'``  *<none>*      clears the cached information
%       =================  ============  ========================================
%
% Output:
%   TorF (boolean) : true if optional functionality is available,
%       false otherwise
%   ver_str (char array) : version number as a string (e.g. ``'3.11.4'``)
%   ver_num (double) : version number as numeric value (e.g. ``3.011004``)
%   rdate (char array) : release date as a string (e.g. ``'29-Feb-2024'``)
%   all_info (struct) : struct with fields:
%
%        - ``av`` -- availability, true if available, false otherwise
%        - ``vstr`` -- version as string
%        - ``vnum`` -- version as numeric value
%        - ``date`` -- release date as string
%
% The optional functionality specified by ``tag`` can be toggled **off**
% or **on** by calling have_feature with a numeric second argument
% ``toggle`` with one of the following values:
%
%   - 0 -- turn **off** the optional functionality
%   - 1 -- turn **on** the optional functionality (if available)
%   - -1 -- **toggle** the on/off state of the optional functionality
%
% Specifying the appropriate string value for the second argument allows
% have_feature to return specific information about the functionality,
% if it is available. This includes the version number as a string or
% numeric value, the release date, or a struct with all of the above.
% For functionality that is not available, all calls with a string-valued
% second argument (except ``'all'``) will return an empty value.
%
% Finally, passing ``'clear_cache'`` as the second argument will cause the
% cached information to be cleared for the specified ``tag`` or, if the first
% argument is ``'all'``, for all optional functionality. When calling with
% ``'clear_cache'`` no return value is defined.
%
% For each valid value of ``tag``, there is a corresponding feature detection
% function named :samp:`have_feature_{<tag>}`, where :samp:`{<tag>}` is the
% ``tag`` value for the feature in question. This makes have_feature modular
% and extensible. Each feature detection function takes no input values, but
% returns three outputs:
%
%   - ``TorF`` -- 1 = feature is available, 0 = feature is not available
%   - ``vstr`` -- version number as a string (e.g. ``'3.11.4'``)
%   - ``rdate`` -- release date as a string (e.g. ``'29-Feb-2024'``)
%
% ``tag`` values for have_feature detection functions included in MP-Test:
%
%   - ``'matlab'`` -- code is running under MATLAB, as opposed to Octave
%   - ``'octave'`` -- code is running under GNU Octave, as opposed to MATLAB
%
% Examples::
%
%   if have_feature('matlab')
%       disp(['Running MATLAB version ', have_feature('matlab', 'vstr')])
%   else
%       disp(['Running Octave version ', have_feature('octave', 'vstr')])
%   end
%
% See also have_feature_matlab, have_feature_octave.

% The following calling syntaxes are also implemented to set and get the
% entire cache struct and are used during testing only.
% ::
%
%     cache = have_feature('all', 'get_cache')
%     have_feature(cache, 'set_cache')

%   MP-Test
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Test.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mptest for more info.

persistent h_f_cache;

action = 'D';                   %% detecting functionality (default)
if nargin > 1
    if ~ischar(rtype) && ~isempty(rtype)
        action = 'T';           %% toggling functionality
        on_off = rtype;
        if on_off < 0                   %% flip the toggle
            TorF = have_feature(tag);
            on_off = ~TorF;
        end
    elseif length(rtype) > 4
        switch lower(rtype)
            case 'get_cache'
                action = 'C';   %% getting cache
                rv = h_f_cache;
            case 'set_cache'
                action = 'C';   %% setting cache
                h_f_cache = tag;
            case 'clear_cache'
                action = 'C';   %% clearing cache
                if strcmpi(tag, 'all')  %% delete all fields
                    h_f_cache = struct();
                else                    %% delete field to force single re-check
                    if isfield(h_f_cache, tag)
                        h_f_cache = rmfield(h_f_cache, tag);
                    end
                end
        end
    end
end

if action == 'T'            %% change availability
    if on_off                   %% turn on if available
        if isfield(h_f_cache, tag)  %% delete field to force single re-check
            h_f_cache = rmfield(h_f_cache, tag);
        end
    else                        %% turn off
        if ~isfield(h_f_cache, tag)     %% not yet been checked
            TorF = have_feature(tag);   %% cache result first
        end
        h_f_cache.(tag).av = 0;         %% then turn off
    end
    TorF = have_feature(tag);           %% return availability
                                        %% (recheck if ON, cached 0 if OFF)
elseif action == 'D'        %% detect availability
    %% info not yet cached?
    if ~isfield(h_f_cache, tag)
        %%-----  determine installation status, version number, etc.  -----
        %% initialize default values
        TorF = 0;
        vstr = '';
        rdate = '';

        %% check for feature
        fcn = ['have_feature_' tag];
        if isempty(which(fcn))
            warning('have_feature: unknown functionality ''%s''', tag);
            vstr = 'unknown';
        else
            [TorF, vstr, rdate] = feval(fcn);
        end

        %% assign values to cache
        h_f_cache.(tag).av   = TorF;
        h_f_cache.(tag).vstr = vstr;
        if isempty(vstr)
            h_f_cache.(tag).vnum = [];
        else
            h_f_cache.(tag).vnum = vstr2num(vstr);
        end
        h_f_cache.(tag).date = rdate;
    end
end

%% extract desired values from cache
if action ~= 'C' || nargout
    if nargin < 2 || action == 'T'
        rv = h_f_cache.(tag).av;
    else
        switch lower(rtype)
            case 'vstr'
                rv = h_f_cache.(tag).vstr;
            case 'vnum'
                rv = h_f_cache.(tag).vnum;
            case 'date'
                rv = h_f_cache.(tag).date;
            case 'all'
                rv = h_f_cache.(tag);
            case {'', 'av'}
                rv = h_f_cache.(tag).av;
        end
    end
end

function num = vstr2num(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
