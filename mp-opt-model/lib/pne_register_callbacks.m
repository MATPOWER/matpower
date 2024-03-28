function reg_cb = pne_register_callbacks(my_cbacks, reg_cb)
% pne_register_callbacks - Register PNE callback functions.
% ::
%
%   REG_CB = PNE_REGISTER_CALLBACKS(MY_CBACKS)
%   REG_CB = PNE_REGISTER_CALLBACKS(MY_CBACKS, REG_CB)
%
%   Registers callback functions to be called by PNES_MASTER.
%
%   Inputs:
%       MY_CBACKS : a callback spec, or cell array of callback specs,
%           where a single callback spec is of the form FCN, {FCN} or
%           {FCN, PRIORITY}, where FCN and PRIORITY are:
%           FCN : function handle to the callback function, see
%               PNE_CALLBACK_DEFAULT for details on calling syntax of
%               callback functions
%           PRIORITY : number that determines order of execution for multiple
%                  callback functions, where higher numbers run first,
%                  default priority is 20, where the standard callbacks
%                  are called with the following priority:
%                       pne_callback_nose       51
%                       pne_callback_target_lam 50
%                       pne_callback_default    0
%       REG_CB : (optional) struct array containing existing registered
%           callback functions
%
%   Outputs:
%       REG_CB : updated struct containing registered callback functions
%
%   User Defined PNES_MASTER Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       PNE_CALLBACK_DEFAULT. These are specified via the 'callbacks' option
%       (e.g. OPT.callbacks) to PNES_MASTER, which takes the same form as
%       MY_CBACKS above.
%
% See also pnes_master, pne_callback_default, pne_callback_nose,
% pne_callback_target_lam.

%   MP-Opt-Model
%   Copyright (c) 2016-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% initialize registered event function struct
if nargin < 2
    reg_cb = [];
end

%% convert input to cell array of cell arrays, if necessary
if ~iscell(my_cbacks)
    my_cbacks = {my_cbacks};
end
for k = 1:length(my_cbacks)
    if ~iscell(my_cbacks{k})
        my_cbacks{k} = {my_cbacks{k}};
    end
end

for k = 1:length(my_cbacks)
    params = my_cbacks{k};      %% callback function parameters
    if length(params) > 1
        priority = params{2};
    else
        priority = 20;
    end
    e = struct( 'fcn', params{1}, 'priority', priority );
    if isempty(reg_cb)
        reg_cb = e;
    else
        reg_cb(end+1) = e;
    end
end

%% sort by descending priority
ncb = length(reg_cb);
p = cell(ncb, 1);
[p{:}] = deal(reg_cb.priority);
[~, i] = sort(cell2mat(p), 'descend');
reg_cb = reg_cb(i);
