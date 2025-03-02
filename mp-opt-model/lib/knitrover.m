function rv = knitrover(varargin)
% knitrover - Prints or returns installed Artelys Knitro version info.
% ::
%
%   knitrover
%   v = knitrover
%   v = knitrover('all')
%
% When called with an output argument and no input argument, knitrover
% returns the current Artelys Knitro version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% knitrover without assigning the return value prints the version and
% release date of the current installation of Artelys Knitro.
%
% See also mpver, knitro_solve.

%   MP-Opt-Model
%   Copyright (c) 2010-2025, Power Systems Engineering Research Center (PSERC)
%   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

kntr = have_feature('knitro', 'all');
if ~kntr.av
    kntr.vstr = '<unknown>';
end

v = struct( 'Name',     'Artelys Knitro', ... 
            'Version',  kntr.vstr, ...
            'Release',  '', ...
            'Date',     kntr.date );
if nargout > 0
    if nargin > 0
        rv = v;
    else
        rv = v.Version;
    end
else
    if kntr.av
        fprintf('%-22s Version %-10s %-11s\n', v.Name, v.Version, v.Date);
    else
        fprintf('%-22s -- not installed --\n', v.Name);
    end
end
