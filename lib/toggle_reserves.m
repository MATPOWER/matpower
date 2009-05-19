function mpc = toggle_reserves(mpc, on_off)
%TOGGLE_RESERVES Enable or disable set of fixed reserves userfcn callbacks.
%   mpc = toggle_reserves(mpc, on_off)
%
%   The 2nd arg (on_off) takes one of two values 'on' or 'off'.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(on_off, 'on')
    %% check for proper reserve inputs
    if ~isfield(mpc, 'reserves') || ~isstruct(mpc.reserves) || ...
            ~isfield(mpc.reserves, 'zones') || ...
            ~isfield(mpc.reserves, 'req') || ...
            ~isfield(mpc.reserves, 'cost')
        error('toggle_reserves: case must contain a ''reserves'' field, a struct defining ''zones'', ''req'' and ''cost''');
    end
    
    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', 'userfcn_reserves_ext2int');
    mpc = add_userfcn(mpc, 'formulation', 'userfcn_reserves_formulation');
    mpc = add_userfcn(mpc, 'int2ext', 'userfcn_reserves_int2ext');
    mpc = add_userfcn(mpc, 'printpf', 'userfcn_reserves_printpf');
    mpc = add_userfcn(mpc, 'savecase', 'userfcn_reserves_savecase');
    
    % mpc = add_userfcn(mpc, 'formulation', 'userfcn_reserves', mpc.reserves);
elseif strcmp(on_off, 'off')
    mpc = remove_userfcn(mpc, 'savecase', 'userfcn_reserves_savecase');
    mpc = remove_userfcn(mpc, 'printpf', 'userfcn_reserves_printpf');
    mpc = remove_userfcn(mpc, 'int2ext', 'userfcn_reserves_int2ext');
    mpc = remove_userfcn(mpc, 'formulation', 'userfcn_reserves_formulation');
    mpc = remove_userfcn(mpc, 'ext2int', 'userfcn_reserves_ext2int');
else
    error('toggle_reserves: 2nd argument must be either ''on'' or ''off''');
end

return;
