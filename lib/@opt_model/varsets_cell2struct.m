function vs = varsets_cell2struct(om, vs)
%VARSETS_CELL2STRUCT  Converts VARSETS from cell array to struct array.
%   VARSETS = OM.VARSETS_CELL2STRUCT(VARSETS)
%
%   Converts VARSETS from a cell array to a struct array, if necessary.
%
%   See also VARSETS_LEN

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% convert varsets from cell to struct array if necessary
if ~isempty(vs) && iscell(vs)
    empty_cells = cell(1, length(vs));
    [empty_cells{:}] = deal({});    %% empty cell arrays
    vs = struct('name', vs, 'idx', empty_cells);
end
