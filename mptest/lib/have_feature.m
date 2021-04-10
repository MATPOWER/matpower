function rv = have_feature(tag, rtype)
%HAVE_FEATURE  Test for optional functionality / version info.
%   TORF = HAVE_FEATURE(TAG)
%   TORF = HAVE_FEATURE(TAG, TOGGLE)
%   VER_STR = HAVE_FEATURE(TAG, 'vstr')
%   VER_NUM = HAVE_FEATURE(TAG, 'vnum')
%   DATE    = HAVE_FEATURE(TAG, 'date')
%   INFO    = HAVE_FEATURE(TAG, 'all')
%   HAVE_FEATURE(TAG, 'clear_cache')
%   HAVE_FEATURE('all', 'clear_cache')
%
%   Returns availability, version and release information for optional
%   functionality. All information is cached, and the cached values
%   returned on subsequent calls. If the functionality exists, an attempt is
%   made to determine the release date and version number. The second
%   argument defines which value is returned, as follows:
%       <none>      1 = optional functionality is available, 0 = not available
%       'vstr'      version number as a string (e.g. '3.11.4')
%       'vnum'      version number as numeric value (e.g. 3.011004)
%       'date'      release date as a string (e.g. '21-Sep-2020')
%       'all'       struct with fields named 'av' (for 'availability'), 'vstr',
%                   'vnum' and 'date', and values corresponding to the above,
%                   respectively.
%
%   For functionality that is not available, all calls with a string-valued
%   second argument will return an empty value.
%
%   Alternatively, the optional functionality specified by TAG can be toggled
%   OFF or ON by calling HAVE_FEATURE with a numeric second argument TOGGLE
%   with one of the following values:
%       0 - turn OFF the optional functionality
%       1 - turn ON the optional functionality (if available)
%      -1 - toggle the ON/OFF state of the optional functionality
%
%   Finally, passing 'clear_cache' as the second argument will cause the
%   cached information to be cleared for the specified TAG or, if the first
%   argument is 'all', for all optional functionality. When calling with
%   'clear_cache' no return value is defined.
%
%   For each valid value of TAG, there is a corresponding feature detection
%   function named HAVE_FEATURE_<TAG>, where <TAG> is the TAG value for the
%   feature in question. This makes HAVE_FEATURE modular and extensible.
%   Each feature detection function takes no input values, but returns
%   three outputs
%       TORF - 1 = feature is available, 0 = feature is not available
%       VSTR - version number as a string (e.g. '3.11.4')
%       RDATE - release date as a string (e.g. '21-Sep-2020')
%
%   TAG values for HAVE_FEATURE detection functions included in MP-Test:
%       matlab      - code is running under MATLAB, as opposed to Octave
%       octave      - code is running under GNU Octave, as opposed to MATLAB
%
%   Examples:
%       if have_feature('matlab')
%           disp(['Running MATLAB version ', have_feature('matlab', 'vstr')])
%       else
%           disp(['Running Octave version ', have_feature('octave', 'vstr')])
%       end
%
%   See also HAVE_FEATURE_MATLAB, HAVE_FEATURE_OCTAVE

%   The following calling syntaxes are also implemented to set and get the
%   entire cache struct and are used during testing only.
%       CACHE = HAVE_FEATURE('all', 'get_cache')
%       HAVE_FEATURE(CACHE, 'set_cache')

%   MP-Test
%   Copyright (c) 2004-2020, Power Systems Engineering Research Center (PSERC)
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
