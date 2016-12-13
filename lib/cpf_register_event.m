function cpf_events = cpf_register_event(cpf_events, name, fcn, tol, locate)
%CPF_REGISTER_EVENT  Register event functions=
%   CPF_EVENTS = CPF_REGISTER_EVENT(CPF_EVENTS, NAME, FCN, TOL, LOCATE)
%
%   Registers a CPF event function to be called by RUNCPF.
%
%   Inputs:
%       CPF_EVENTS : struct containing info about registered CPF event fcns
%       NAME : string containing event name
%       FCN : string containing name of event function, returning numerical
%             scalar or vector value that changes sign at location of the event
%       TOL : scalar or vector of same dimension as event function return value
%             of tolerance for detecting the event, i.e. abs(val) <= tol
%       LOCATE : flag indicating whether the event requests a rollback step
%                to locate the event function zero
%
%   Outputs:
%       CPF_EVENTS : updated struct containing info about registered CPF event fcns

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Shrirang Abhyankar, Argonne National Laboratory
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% the event function data to be registered
e = struct( ...
        'name', name, ...
        'fcn', fcn, ...
        'tol', tol, ...
        'locate', locate ...
    );

%% convert function names to function handles, as necessary
if ~isa(e.fcn, 'function_handle')
    e.fcn = str2func(e.fcn);
end

%% register to list of event functions
if isempty(cpf_events)
    cpf_events = e;
else
    nef = length(cpf_events);
    for k = 1:nef
        if strcmp(cpf_events(k).name, name)
            error('cpf_register_event: duplicate event name: ''%s''', name);
        end
    end
    cpf_events(nef+1) = e;
end
