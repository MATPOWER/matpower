function sd = loadstoragedata(sd_table, mpc)
%LOADSTORAGEDATA   Load StorageData table into StorageData struct.
%
%   SD = LOADSTORAGEDATA(SD_TABLE)
%   SD = LOADSTORAGEDATA(SD_TABLE, GEN)
%   SD = LOADSTORAGEDATA(SD_TABLE, MPC)
%
%   Loads data from a StorageData table struct, or from a M-file or MAT-file
%   that returns such a struct, and converts it to a StorageData struct.
%   StorageData contains all parameters required by MOST for storage
%   units that is not included in MPC or xGenData.
%
%   If the first argument is the name of an M-file the optional second
%   argument will be passed to it. This can be useful for cases where
%   the user wishes to use data from GEN or GENCOST to set values in the
%   StorageData.
%
%   Inputs:
%       SD_TABLE : a StorageData table struct or the name of an M-file
%                   or MAT-file that returns one, with the following fields
%           .colnames : N dimensional cell array of names corresponding
%               to the columns of the 'data' field. Valid column names
%               are the same as the output field names. All columns are
%               optional and the correpsonding defaults are listed
%               next the output field names below.
%           .data : (NS x N) matrix of data
%           .MinStorageLevel                 (optional)     (scalar)
%           .MaxStorageLevel                 (optional)     (scalar)
%           .OutEff                          (optional)     (scalar)
%           .InEff                           (optional)     (scalar)
%           .LossFactor                      (optional)     (scalar)
%           .rho                             (optional)     (scalar)
%           (values in any of the scalar fields listed above are
%            overridden by any corresponding values in the 'data' table)
%       GEN : (optional) standard GEN matrix for generators corresponding
%             to storage units in SD_TABLE
%       MPC : (optional) MATPOWER case struct containing GEN and GENCOST
%             matrices for generators corresponding to storage units in
%             SD_TABLE
%
%   Output:
%       SD : a StorageData struct, contains the following fields,
%               all of which are (NS x 1) vectors unless otherwise
%               indicated. Defaults are as provided by MD_INIT.
%           .UnitIdx                            (ns x 1)
%           .ExpectedTerminalStorageAim         (ns x 1)
%           .ExpectedTerminalStorageMin         (ns x 1)
%           .ExpectedTerminalStorageMax         (ns x 1)
%           .InitialStorage                     (ns x 1)
%           .InitialStorageLowerBound           (ns x 1)
%           .InitialStorageUpperBound           (ns x 1)
%           .InitialStorageCost                 (ns x 1)
%           .TerminalStoragePrice               (ns x 1)
%           .TerminalChargingPrice0             (ns x 1)
%           .TerminalDischargingPrice0          (ns x 1)
%           .TerminalChargingPriceK             (ns x 1)
%           .TerminalDischargingPriceK          (ns x 1)
%           .MinStorageLevel                    (ns x 1) or (1 x 1)
%           .MaxStorageLevel                    (ns x 1) or (1 x 1)
%           .OutEff                             (ns x 1) or (1 x 1)
%           .InEff                              (ns x 1) or (1 x 1)
%           .LossFactor                         (ns x 1) or (1 x 1)
%           .rho                                (ns x 1) or (1 x 1)

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
if ischar(sd_table)
    infile = sprintf(' in file: ''%s''', sd_table);
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
ns0 = size(gen, 1);      %% number of storage units in gen matrix

%% load the table
fields = {'colnames', 'data'};
sdt = loadgenericdata(sd_table, 'struct', fields, 'sd_table', args);
if isempty(sdt)
    sdt.colnames = {};
    sdt.data = [];
end
[ns, nc] = size(sdt.data);

%% check consistency of dimensions of data table, colnames and gen
if ns0 > 0 && ns > 0 && ns ~= ns0
    error('loadstoragedata: # of rows in ''data'' table (%d) do not equal rows in GEN (%d)%s', ...
        ns, ns0, infile);
end
if nc ~= length(sdt.colnames)
    error('loadstoragedata: # of columns in ''data'' table (%d) do match entries in ''colnames'' (%d)%s', ...
        nc, length(sdt.colnames), infile);
end

%% initialize storage data with default values from md_init's Storage field
sd = getfield(md_init(), 'Storage');

%% set individual fields
f = {   'MinStorageLevel', ...
        'MaxStorageLevel', ...
        'OutEff', ...
        'InEff', ...
        'LossFactor', ...
        'rho'
    };
for k = 1:length(f)
    if isfield(sd_table, f{k});
        sd.(f{k}) = sd_table.(f{k});
    end
end

%% create cell array with valid fields for checking
valid_fields = {
        'UnitIdx', ...
        'ExpectedTerminalStorageAim', ...
        'ExpectedTerminalStorageMin', ...
        'ExpectedTerminalStorageMax', ...
        'InitialStorage', ...
        'InitialStorageLowerBound', ...
        'InitialStorageUpperBound', ...
        'InitialStorageCost', ...
        'TerminalStoragePrice', ...
        'TerminalChargingPrice0', ...
        'TerminalDischargingPrice0', ...
        'TerminalChargingPriceK', ...
        'TerminalDischargingPriceK', ...
        'MinStorageLevel', ...
        'MaxStorageLevel', ...
        'OutEff', ...
        'InEff', ...
        'LossFactor', ...
        'rho'
    };

%% copy data from table
for j = 1:nc
    %% check for valid colname
    if ~any(strcmp(sdt.colnames{j}, valid_fields))
        error('loadstoragedata: ''%s'' is not a valid StorageData field name%s', sdt.colnames{j}, infile);
    end
    %% copy data
    sd.(sdt.colnames{j}) = sdt.data(:, j);
end
