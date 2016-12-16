function xgd_table = ex_xgd(mpc)
%EX_XGD  Example xGenData table for stochastic OPF.

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
    1e-8    250     2e-8    250     1e-9    1e-9;
    1e-8    250     2e-8    250     1e-9    1e-9;
    1e-8    600     2e-8    600     1e-9    1e-9;
    1e-8    800     2e-8    800     1e-9    1e-9;
];
