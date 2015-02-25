function val = get(om, varargin)
%GET  Returns the value of a field.
%   VAL = GET(OM, FIELD1, FIELD2, ...)
%
%   Example:
%       var_order = get(om, 'var', 'order');
%
%   See also OPT_MODEL.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008-2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://matpower.org/ for more info.

val = om;
for k = 1:length(varargin)
    if ischar(varargin{k})
        val = val.(varargin{k});
    else
        val = val(varargin{k});
    end
end
