function varargout = qps_matpower(varargin)
% qps_matpower - Deprecated, please use qps_master instead.

%   MATPOWER
%   Copyright (c) 2010-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

[varargout{1:nargout}] = qps_master(varargin{:});
