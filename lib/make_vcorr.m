function V_corr = make_vcorr(DD,pv,nb,nl,f,Zb)
% make_vcorr - Voltage Correction used in distribution power flow.
% ::
%
%   V_corr = make_vcorr(DD,pv,nb,nl,f,Zb)
%
%   Calculates voltage corrections with current generators placed at PV
%   buses. Their currents are calculated with the voltage difference at PV
%   buses break points and loop impedances. The slack bus voltage is set to
%   zero. Details can be seen in
%   D. Rajicic, R. Ackovski and R. Taleski, "Voltage correction power flow,"
%   IEEE Transactions on Power Delivery, vol. 9, no. 2, pp. 1056-1062, Apr 1994.
%   https://doi.org/10.1109/61.296308
%
% See also radial_pf.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Mirko Todorovski
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

V_corr = zeros(nb,1);
I = zeros(nb,1);
I(pv) = DD;
% backward sweep
for k = nl:-1:2
    i = f(k);
    I(i) = I(i) + I(k);
end
% forward sweep
for k = 2:nl
    i = f(k);
    V_corr(k) = V_corr(i) - Zb(k) * I(k);
end
