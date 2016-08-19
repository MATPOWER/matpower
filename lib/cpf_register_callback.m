function cpf_callbacks = cpf_register_callback(cpf_callbacks, name, fcn)
%CPF_REGISTER_CALLBACK  Register CPF callback functions
%
%   CPF_CALLBACKS = CPF_REGISTER_CALLBACK(CPF_CALLBACKS, NAME, FCN)
%   
%   Inputs:
%       CPF_CALLBACKS : struct containing info about registered CPF
%                       callback fcns
%       NAME : string containing callback identifier or name
%       FCN : string containing name of callback function
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

cb = struct( ...
        'name', name, ...
        'fcn', fcn ...
    );
if ~isa(cb.fcn, 'function_handle')
    cb.fcn = str2func(cb.fcn);
end
if isempty(cpf_callbacks)
    cpf_callbacks = cb;
else
    cpf_callbacks(end+1) = cb;
end
