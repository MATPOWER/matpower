function [TorF, vstr, rdate] = have_feature_ipopt_auxdata()
% have_feature_ipopt_auxdata - Detect availability/version info for :func:`ipopt_auxdata`.
%
% Private feature detection function implementing ``'ipopt_auxdata'`` tag
% for have_feature to detect support for :func:`ipopt_auxdata`, required
% by IPOPT 3.11 and later.
%
% See also have_feature, nlps_master, qps_master, ipopt, ipopt_auxdata.

%   MP-Opt-Model
%   Copyright (c) 2004-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

vstr = '';
rdate = '';
if have_feature('ipopt')
    vn = have_feature('ipopt', 'vnum');
    if ~isempty(vn)
        if vn >= 3.011  %% have_feature('ipopt') already checked
            TorF = 1;   %% for existence of 'ipopt_auxdata'
        else
            TorF = 0;   %% don't use it, even if it exists
        end
    %% no version info, decide based on presence
    %% or absence of 'ipopt_auxdata'
    elseif exist('ipopt_auxdata', 'file')
        TorF = 1;
    else
        TorF = 0;
    end
else
    TorF = 0;
end
