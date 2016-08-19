function cpf_callbacks = cpf_register_callback(cpf_callbacks, fcn, priority)
%CPF_REGISTER_CALLBACK  Register CPF callback functions
%
%   CPF_CALLBACKS = CPF_REGISTER_CALLBACK(CPF_CALLBACKS, FCN, PRIORITY)
%   
%   Inputs:
%       CPF_CALLBACKS : struct containing info about registered CPF
%                       callback fcns
%       FCN : string containing name of callback function
%       PRIORITY : number that determines order of execution for multiple
%                  callback functions, where higher numbers run first,
%                  default priority is 20, where the standard callbacks
%                  are called with the following priority:
%                       cpf_nose_event_cb       51
%                       cpf_target_lam_event_cb 50
%                       cpf_qlim_event_cb       41
%                       cpf_plim_event_cb       40
%                       cpf_default_callback    0
%
%   Outputs:
%       CPF_CALLBACKS : updated struct containing info about registered
%                       CPF callback fcns

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default priority
if nargin < 3 || isempty(priority)
    priority = 20;
end

cb = struct( ...
        'fcn', fcn, ...
        'priority', priority ...
    );
if ~isa(cb.fcn, 'function_handle')
    cb.fcn = str2func(cb.fcn);
end

%% add it to the list
if isempty(cpf_callbacks)
    cpf_callbacks = cb;         %% first one
else
    ncb = length(cpf_callbacks) + 1;
    cpf_callbacks(ncb) = cb;    %% append
    
    %% sort by descending priority
    p = cell(ncb, 1);
    [p{:}] = deal(cpf_callbacks.priority);
    [junk, i] = sort(cell2mat(p), 'descend');
    cpf_callbacks = cpf_callbacks(i);
end
