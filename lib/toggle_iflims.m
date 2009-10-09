function mpc = toggle_iflims(mpc, on_off)
%TOGGLE_IFLIMS Enable or disable set of interface flow limit userfcn callbacks.
%   mpc = toggle_iflims(mpc, on_off)
%
%   The 2nd arg (on_off) takes one of two values 'on' or 'off'.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if strcmp(on_off, 'on')
    %% check for proper reserve inputs
    if ~isfield(mpc, 'if') || ~isstruct(mpc.if) || ...
            ~isfield(mpc.if, 'map') || ...
            ~isfield(mpc.if, 'lims')
        error('toggle_iflims: case must contain an ''if'' field, a struct defining ''map'' and ''lims''');
    end
    
    %% add callback functions
    %% note: assumes all necessary data included in 1st arg (mpc, om, results)
    %%       so, no additional explicit args are needed
    mpc = add_userfcn(mpc, 'ext2int', 'userfcn_iflims_ext2int');
    mpc = add_userfcn(mpc, 'formulation', 'userfcn_iflims_formulation');
    mpc = add_userfcn(mpc, 'int2ext', 'userfcn_iflims_int2ext');
    mpc = add_userfcn(mpc, 'printpf', 'userfcn_iflims_printpf');
    mpc = add_userfcn(mpc, 'savecase', 'userfcn_iflims_savecase');
elseif strcmp(on_off, 'off')
    mpc = remove_userfcn(mpc, 'savecase', 'userfcn_iflims_savecase');
    mpc = remove_userfcn(mpc, 'printpf', 'userfcn_iflims_printpf');
    mpc = remove_userfcn(mpc, 'int2ext', 'userfcn_iflims_int2ext');
    mpc = remove_userfcn(mpc, 'formulation', 'userfcn_iflims_formulation');
    mpc = remove_userfcn(mpc, 'ext2int', 'userfcn_iflims_ext2int');
else
    error('toggle_iflims: 2nd argument must be either ''on'' or ''off''');
end

return;
