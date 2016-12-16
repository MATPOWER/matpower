function [idx, new_mpc, new_xgd, new_sd] = addstorage(storage, mpc, xgd, sd)
%ADDWIND   Adds storage units and corresponding xGen/StorageData to existing data.
%
%   [IDX, NEW_MPC] = ADDSTORAGE(STORAGE, MPC)
%   [IDX, NEW_MPC, NEW_XGD, NEW_SD] = ADDSTORAGE(STORAGE, MPC)
%   [IDX, NEW_MPC, NEW_XGD, NEW_SD] = ADDSTORAGE(STORAGE, MPC, XGD)
%   [IDX, NEW_MPC, NEW_XGD, NEW_SD] = ADDSTORAGE(STORAGE, MPC, XGD, SD)
%
%   Given a StorageUnitData structure, or the name of a file containing such
%   a structure, this function adds the specified storage generators to an
%   existing MATPOWER case, xGenData struct and StorageData struct.
%
%   Inputs:
%       STORAGE : a StorageUnitData struct or the name of an M-file
%                   or MAT-file that returns one, with the following fields
%           .gen     : rows to be appended to the GEN matrix from MPC
%           .gencost : (optional) rows to be added to the GENCOST matrix
%                      from MPC, default is zero cost
%           .xgd_table : xGenData table struct or filename providing data for
%                        the storage units being added. See LOADSTORAGEDATA
%                        for more information on the xGenData table format.
%           .sd_table : StorageData table struct or filename providing data
%                       for the storage units being added. See LOADSTORAGEDATA
%                       for more information on the StorageData table format.
%       MPC : MATPOWER case struct to which storage generators will be added
%       XGD : (optional) xGenData struct corresponding to the generators
%             already in MPC, to which the new xGenData for the storage
%             units will be added.
%       SD : (optional) StorageData struct corresponding to the generators
%            already in MPC, to which the new StorageData for the storage
%            units will be added.
%
%   Output:
%       IDX : Generator indices of newly added storage units.
%       NEW_MPC : MPC with storage units appended to MPC.GEN and MPC.GENCOST
%                 MPC.GENFUEL (= 'ess').
%       NEW_XGD : XGD with xGenData for new storage units appended.
%       NEW_SD  : SD with StorageData for new storage units appended.
%
%   See also LOADSTORAGEDATA, LOADXGENDATA.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% define fuel type for storage
STORAGE_FUEL = 'ess';

%% input arg handling
if nargin < 4
    sd = [];
    if nargin < 3
        xgd = [];
    end
end
if ischar(storage)
    infile = sprintf(' in file: ''%s''', storage);
    storage = loadgenericdata(storage, 'struct', 'gen', 'storage');
else
    infile = '';
end

%% add to MPC
ns = size(storage.gen, 1);          %% number of storage units being added
if isfield(storage, 'gencost')
    storage_gencost = storage.gencost;
else        %% use zero cost by default
    storage_gencost = repmat([POLYNOMIAL 0 0 2 0 0], ns, 1);
end
[new_mpc, idx] = addgen2mpc(mpc, storage.gen, storage_gencost, STORAGE_FUEL);

%% handle xGenData and StorageData
if nargout > 2      %% output NEW_XGD, NEW_SD requested
    %% xGenData
    if isfield(storage, 'xgd_table')
        storage_xgd = loadxgendata(storage.xgd_table, storage.gen);
    else
        error('addstorage: missing XGD_TABLE field in STORAGE');
    end

    if isempty(xgd)     %% no input XGD provided
        new_xgd = storage_xgd;
    else                %% input XGD provided
        %% append rows of every field in xgd
        new_xgd = xgd;
        fields = fieldnames(xgd);
        for f = 1:length(fields)
            ff = fields{f};
            %% dims of wind_xgd fields already checked by loadxgendata
            if size(xgd.(ff), 1) ~= size(mpc.gen, 1)
                error('addstorage: # of rows in XGD.%s (%d) does not match MPC.GEN (%d)', ...
                    ff, size(xgd.(ff), 1), size(mpc.gen, 1));
            end
            new_xgd.(ff) = [xgd.(ff); storage_xgd.(ff)];
        end
    end

    %% StorageData
    if isfield(storage, 'sd_table')
        storage_sd = loadstoragedata(storage.sd_table, storage.gen);
    else
        error('addstorage: missing SD_TABLE field in STORAGE');
    end

    storage_sd.UnitIdx = idx;       %% add UnitIdx field

    if isempty(sd)     %% no input SD provided
        new_sd = storage_sd;
    else                %% input SD provided
        %% find number of storage units in original mpc
        ns0 = length(find(strcmp(mpc.genfuel, STORAGE_FUEL)));
        new_sd = sd;
        fields = fieldnames(sd);
        for f = 1:length(fields)    %% append rows of every field in sd
            ff = fields{f};
            add_ns  = size(storage_sd.(ff), 1);
            orig_ns = size(        sd.(ff), 1);
            if orig_ns ~= 1 || add_ns ~= 1 || sd.(ff) ~= storage_sd.(ff)
                %% at least one of them is not 1-d or else they don't match
                if orig_ns == 1     %% original is 1-d
                    %% expand it
                    sd.(ff) = ones(ns0, 1) * sd.(ff);
                end
                if add_ns == 1      %% new is 1-d
                    %% expand it
                    storage_sd.(ff) = ones(ns, 1) * storage_sd.(ff);
                end
                %% concatenate them
                new_sd.(ff) = [sd.(ff); storage_sd.(ff)];
            % else (both are 1-d and they match)
            %   do nothing, just use old sd.(ff)
            end
        end
    end
end
