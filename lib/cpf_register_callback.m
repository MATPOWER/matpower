function cpf_callbacks = cpf_register_callback(cpf_callbacks, fcn, priority, args)
%CPF_REGISTER_CALLBACK  Register CPF callback functions
%   CPF_CALLBACKS = CPF_REGISTER_CALLBACK(CPF_CALLBACKS, FCN, PRIORITY)
%
%   Registers a CPF callback function to be called by RUNCPF.
%
%   Inputs:
%       CPF_CALLBACKS : struct containing info about registered CPF
%                       callback fcns
%       FCN : string containing name of callback function
%       PRIORITY : number that determines order of execution for multiple
%                  callback functions, where higher numbers run first,
%                  default priority is 20, where the standard callbacks
%                  are called with the following priority:
%                       cpf_flim_event_cb       53
%                       cpf_vlim_event_cb       52
%                       cpf_nose_event_cb       51
%                       cpf_target_lam_event_cb 50
%                       cpf_qlim_event_cb       41
%                       cpf_plim_event_cb       40
%                       cpf_default_callback    0
%^      ARGS : arguments to be passed to the callback each time it is invoked
%
%   Outputs:
%       CPF_CALLBACKS : updated struct containing info about registered
%                       CPF callback fcns
%
%   User Defined CPF Callback Functions:
%       The user can define their own callback functions which take
%       the same form and are called in the same contexts as
%       CPF_DEFAULT_CALLBACK. These are specified via the MATPOWER
%       option 'cpf.user_callback'. This option can be a string containing
%       the name of the callback function, or a struct with the following
%       fields, where all but the first are optional:
%           'fcn'       - string with name of callback function
%           'priority'  - numerical value specifying callback priority
%                (default = 20, see CPF_REGISTER_CALLBACK for details)
%           'args'      - arbitrary value (any type) passed to the callback
%                         as CB_ARGS each time it is invoked
%       Multiple user callbacks can be registered by assigning a cell array
%       of such strings and/or structs to the 'cpf.user_callback' option.
%
%   See also RUNCPF, CPF_DEFAULT_CALLBACK.

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default priority
if nargin < 4
    args = [];
    if nargin < 3
        priority = [];
    end
end
if isempty(priority)
    priority = 20;
end

cb = struct( ...
        'fcn', fcn, ...
        'priority', priority, ...
        'args', args ...
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
