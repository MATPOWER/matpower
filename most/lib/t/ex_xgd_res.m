function xgd_table = ex_xgd_res(mpc)
%EX_XGD_RES  Example xGenData table for stochastic OPF w/reserve costs.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% initial xGenData
xgd_table.colnames = {
    'PositiveActiveReservePrice', ...
            'PositiveActiveReserveQuantity', ...
                    'NegativeActiveReservePrice', ...
                            'NegativeActiveReserveQuantity', ...
                                    'PositiveActiveDeltaPrice', ...
                                            'NegativeActiveDeltaPrice', ...
};
xgd_table.data = [
    5       250     10      0       1e-9    1e-9;
    1e-8    100     2e-8    0       1e-9    1e-9;
    1.5     600     3       0       1e-9    1e-9;
    1e-8    800     2e-8    0       1e-9    1e-9;
];
