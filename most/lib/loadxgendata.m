function xgd = loadxgendata(xgd_table, mpc)
%LOADXGENDATA   Load xGenData table into xGenData struct.
%
%   XGD = LOADXGENDATA(XGD_TABLE)
%   XGD = LOADXGENDATA(XGD_TABLE, GEN)
%   XGD = LOADXGENDATA(XGD_TABLE, MPC)
%
%   Loads data from an xGenData table struct, or from a function or MAT-file
%   that returns such a struct, and converts it to an xGenData struct.
%   xGenData contains all of the per-generator data required by MOST
%   that is not included in MPC, including reserve offer data.
%
%   If the first argument is the name of a function the optional second
%   argument will be passed to it. This can be useful for cases where
%   the user wishes to use data from GEN or GENCOST to set values in the
%   xGenData.
%
%   Inputs:
%       XGD_TABLE : an xGenData table struct or the name of a function
%                   or MAT-file that returns one, with 2 fields
%           .colnames : N dimensional cell array of names corresponding
%               to the columns of the 'data' field. Valid column names
%               are the same as the output field names. All columns are
%               optional and the corresponding defaults are listed
%               next to the output field names below.
%           .data : (NG x N) matrix of data
%           If XGD_TABLE is empty, an xGenData struct will be created
%           with default values (requires GEN or MPC).
%       GEN : (optional) standard GEN matrix for generators corresponding
%             to XGD_TABLE
%       MPC : (optional) MATPOWER case struct containing GEN and GENCOST
%             matrices for generators corresponding to XGD_TABLE
%
%   Output:
%       XGD :   an xGenData struct, contains the following fields,
%               all of which are (NG x 1) vectors
%           .CommitSched                        (default = C*)
%               (if MPC or GEN not provided, default = 1)
%           .InitialPg                          (default = GEN(:, PG))
%               (only optional on input if MPC or GEN are provided)
%           .TerminalPg (only included if provided in input)
%           .RampWearCostCoeff                  (default = 0)
%           .PositiveActiveReservePrice         (default = 0)
%           .PositiveActiveReserveQuantity      (default = R*)
%           .NegativeActiveReservePrice         (default = 0)
%           .NegativeActiveReserveQuantity      (default = R*)
%           .PositiveActiveDeltaPrice           (default = 0)
%           .NegativeActiveDeltaPrice           (default = 0)
%           .PositiveLoadFollowReservePrice     (default = 0)
%           .PositiveLoadFollowReserveQuantity  (default = R*)
%           .NegativeLoadFollowReservePrice     (default = 0)
%           .NegativeLoadFollowReserveQuantity  (default = R*)
%           .CommitKey                          (default = not present)
%           .InitialState **                    (default = +/-Inf, based on C*)
%           .MinUp **                           (default = 1)
%           .MinDown **                         (default = 1)
%               (potential future additions)
%                   FuelType (would need to define numeric constants for each)
%
%       * If GEN or MPC are provided then C = GEN(:, GEN_STATUS) and
%         R = 2 * (GEN(:, PMAX) - MIN(0, GEN(:, PMIN))), otherwise
%         C = 1 and R = Inf.
%       ** Requires CommitKey be present and non-empty.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% input arg handling
if ischar(xgd_table)
    infile = sprintf(' in file: ''%s''', xgd_table);
else
    infile = '';
end
if nargin < 2
    args = [];
    gen = [];
else
    args = mpc;
    if isstruct(mpc)
        gen = mpc.gen;
    else
        gen = mpc;
    end
end
ng0 = size(gen, 1);      %% number of gens in gen matrix

%% load the table
fields = {'colnames', 'data'};
xgdt = loadgenericdata(xgd_table, 'struct', fields, 'xgd_table', args);
if isempty(xgdt)
    xgdt.colnames = {};
    xgdt.data = [];
end
[ng, nc] = size(xgdt.data);

%% check consistency of dimensions of data table, colnames and gen
if ng0 > 0 && ng > 0 && ng ~= ng0
    error('loadxgendata: # of rows in ''data'' table (%d) do not equal rows in GEN (%d)%s', ...
        ng, ng0, infile);
end
if nc ~= length(xgdt.colnames)
    error('loadxgendata: # of columns in ''data'' table (%d) do match entries in ''colnames'' (%d)%s', ...
        nc, length(xgdt.colnames), infile);
end

%% check for UC
if any(strcmp('CommitKey', xgdt.colnames))
    UC = 1;
else
    UC = 0;
end

%% initialize xGenData output struct with defaults
if ng0 > 0
    C = gen(:, GEN_STATUS);
    Pg = gen(:, PG);
    R = 2 * (gen(:, PMAX) - min(0, gen(:, PMIN)));
    Z = zeros(ng0, 1);
    W = ones(ng0, 1);
else
    C = ones(ng, 1);
    Pg = Z;
    R = Inf(ng, 1);
    Z = zeros(ng, 1);
    W = ones(ng, 1);
end
xgd = struct( ...
    'CommitSched',                          C, ...
    'InitialPg',                            Pg, ...
    'RampWearCostCoeff',                    Z, ...
    'PositiveActiveReservePrice',           Z, ...
    'PositiveActiveReserveQuantity',        R, ...
    'NegativeActiveReservePrice',           Z, ...
    'NegativeActiveReserveQuantity',        R, ...
    'PositiveActiveDeltaPrice',             Z, ...
    'NegativeActiveDeltaPrice',             Z, ...
    'PositiveLoadFollowReservePrice',       Z, ...
    'PositiveLoadFollowReserveQuantity',    R, ...
    'NegativeLoadFollowReservePrice',       Z, ...
    'NegativeLoadFollowReserveQuantity',    R ...
);
if UC
    xgd.CommitKey       = C;
    xgd.InitialState    = C;
    xgd.InitialState(C > 0)  =  Inf;
    xgd.InitialState(C <= 0) = -Inf;
    xgd.MinUp           = W;
    xgd.MinDown         = W;
end

%% create cell array with valid fields for checking
valid_fields = fieldnames(xgd);
valid_fields{end+1} = 'TerminalPg';
valid_fields{end+1} = 'CommitKey';
valid_fields{end+1} = 'InitialState';
valid_fields{end+1} = 'MinUp';
valid_fields{end+1} = 'MinDown';

%% copy data from table
for j = 1:nc
    %% check for valid colname
    if ~any(strcmp(xgdt.colnames{j}, valid_fields))
        error('loadxgendata: ''%s'' is not a valid xGenData field name%s', xgdt.colnames{j}, infile);
    end
    %% copy data
    xgd.(xgdt.colnames{j}) = xgdt.data(:, j);
end
