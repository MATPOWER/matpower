function [idx, new_mpc, new_xgd] = addwind(wind, mpc, xgd)
%ADDWIND   Adds wind generators and corresponding xGenData to existing data.
%
%   [IDX, NEW_MPC] = ADDWIND(WIND, MPC)
%   [IDX, NEW_MPC, NEW_XGD] = ADDWIND(WIND, MPC)
%   [IDX, NEW_MPC, NEW_XGD] = ADDWIND(WIND, MPC, XGD)
%
%   Given a WindUnitData structure, or the name of a file containing such
%   a structure, this function adds the specified wind generators to an
%   existing MATPOWER case and xGenData struct.
%
%   Inputs:
%       WIND : a WindUnitData struct or the name of an M-file
%                   or MAT-file that returns one, with the following fields
%           .gen     : rows to be appended to the GEN matrix from MPC
%           .gencost : (optional) rows to be added to the GENCOST matrix
%                      from MPC, default is zero cost
%           .xgd_table : xGenData table struct or filename providing data for
%                        the wind units being added. See LOADXGENDATA for
%                        more information on the xGenData table format.
%       MPC : MATPOWER case struct to which wind generators will be added
%       XGD : (optional) xGenData struct corresponding to the generators
%             already in MPC, to which the new xGenData for the wind units
%             will be added.
%
%   Output:
%       IDX : Generator indices of newly added wind units.
%       NEW_MPC : MPC with wind units appended to MPC.GEN and MPC.GENCOST
%                 MPC.GENFUEL (= 'wind').
%       NEW_XGD : XGD with xGenData for new wind units appended.
%
%   See also LOADXGENDATA.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% define fuel type for wind
WIND_FUEL = 'wind';

%% input arg handling
if nargin < 3
    xgd = [];
end
if ischar(wind)
    infile = sprintf(' in file: ''%s''', wind);
    wind = loadgenericdata(wind, 'struct', 'gen', 'wind');
else
    infile = '';
end

%% add to MPC
nw = size(wind.gen, 1);
if isfield(wind, 'gencost')
    wind_gencost = wind.gencost;
else        %% use zero cost by default
    wind_gencost = repmat([POLYNOMIAL 0 0 2 0 0], nw, 1);
end
[new_mpc, idx] = addgen2mpc(mpc, wind.gen, wind_gencost, WIND_FUEL);

%% handle xGenData
if nargout > 2      %% output NEW_XGD requested
    if isfield(wind, 'xgd_table')
        wind_xgd = loadxgendata(wind.xgd_table, wind.gen);
    else
        error('addwind: missing XGD_TABLE field in WIND');
    end

    if isempty(xgd)     %% no input XGD provided
        new_xgd = wind_xgd;
    else                %% input XGD provided
        new_xgd = xgd;
        fields = fieldnames(xgd);
        for f = 1:length(fields)    %% append rows of every field in xgd
            ff = fields{f};
            %% dims of wind_xgd fields already checked by loadxgendata
            if size(xgd.(ff), 1) ~= size(mpc.gen, 1)
                error('addwind: # of rows in XGD.%s (%d) does not match MPC.GEN (%d)', ...
                    ff, size(xgd.(ff), 1), size(mpc.gen, 1));
            end
            new_xgd.(ff) = [xgd.(ff); wind_xgd.(ff)];
        end
    end
end
