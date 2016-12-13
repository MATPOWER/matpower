function Sd = makeSdzip(baseMVA, bus, mpopt)
%MAKESDZIP   Builds vectors of nominal complex bus power demands for ZIP loads.
%   SD = MAKESDZIP(BASEMVA, BUS, MPOPT) returns a struct with three fields,
%   each an nb x 1 vectors. The fields 'z', 'i' and 'p' correspond to the
%   nominal p.u. complex power (at 1 p.u. voltage magnitude) of the constant
%   impedance, constant current, and constant power portions, respectively of
%   the ZIP load model.
%
%   Example:
%       Sd = makeSdzip(baseMVA, bus, mpopt);

%   MATPOWER
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Shrirang Abhyankar
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

if nargin < 3
    mpopt = [];
end
if ~isempty(mpopt) && ~isempty(mpopt.exp.sys_wide_zip_loads.pw)
    if any(size(mpopt.exp.sys_wide_zip_loads.pw) ~= [1 3])
        error('makeSdzip: ''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
    end
    if abs(sum(mpopt.exp.sys_wide_zip_loads.pw) - 1) > eps
        error('makeSdzip: elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
    end
    pw = mpopt.exp.sys_wide_zip_loads.pw;
else
    pw = [1 0 0];
end
if ~isempty(mpopt) && ~isempty(mpopt.exp.sys_wide_zip_loads.qw)
    if any(size(mpopt.exp.sys_wide_zip_loads.qw) ~= [1 3])
        error('makeSdzip: ''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
    end
    if abs(sum(mpopt.exp.sys_wide_zip_loads.qw) - 1) > eps
        error('makeSdzip: elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
    end
    qw = mpopt.exp.sys_wide_zip_loads.qw;
else
    qw = pw;
end

Sd.z = (bus(:, PD) * pw(3)  + 1j * bus(:, QD) * qw(3)) / baseMVA;
Sd.i = (bus(:, PD) * pw(2)  + 1j * bus(:, QD) * qw(2)) / baseMVA;
Sd.p = (bus(:, PD) * pw(1)  + 1j * bus(:, QD) * qw(1)) / baseMVA;
