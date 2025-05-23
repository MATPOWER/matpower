function rv = highsver(varargin)
% highsver - Prints or returns installed HiGHS version info.
% ::
%
%   highsver
%   v = highsver
%   v = highsver('all')
%
% When called with an output argument and no input argument, highsver
% returns the current HiGHS version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% highsver without assigning the return value prints the version and
% release date of the current installation of HiGHS.
%
% See also mpver, callhighs.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

h = have_feature('highs', 'all');
if ~h.av
    h.vstr = '<unknown>';
end

v = struct( 'Name',     'HiGHS', ...
            'Version',  h.vstr, ...
            'Release',  '', ...
            'Date',     h.date );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    if h.av
        fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
    else
        fprintf('%-22s -- not installed --\n', v.Name);
    end
end
