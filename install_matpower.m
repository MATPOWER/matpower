function succ = install_matpower(modify, save_it, verbose, rm_oldpaths)
% install_matpower - Assist user in setting path to install |MATPOWER|.
% ::
%
%   install_matpower
%   install_matpower(modify)
%   install_matpower(modify, save_it)
%   install_matpower(modify, save_it, verbose)
%   install_matpower(modify, save_it, verbose, rm_oldpaths)
%   success = install_matpower(...)
%
% Assists the user in setting up the proper |MATLAB|/Octave path to
% be able to use |MATPOWER| and run its tests. With no input arguments
% it prompts interactively to determine how to handle the paths.
%
% .. note::
%
%   This function is generally **not** in your |MATLAB|/Octave path, unless
%   you change your current working directory to the |MATPOWER| install
%   directory where it is located.
%
% There are two main approaches for installing |MATPOWER|.
%
%   1. If you have a single version of |MATPOWER|, select the options to
%      modify and save the path (interactive option 3). This will add
%      |MATPOWER| to your default |MATLAB|/Octave path for all future
%      sessions.
%
%   2. If you have multiple versions of |MATPOWER|, select the options to
%      not modify the path (interactive option 1), but to save the
%      :func:`addpath` commands to a file. Then execute the saved file to
%      use this version of |MATPOWER|.
%
% *All inputs and outputs are optional.*
%
% Inputs:
%   modify (boolean) : select how to set path
%
%       - 0 *(default)* - generate relevant :func:`addpath` commands, but
%         don't execute them; |MATPOWER| is not installed
%       - 1 - modify the path by executing the relevant :func:`addpath`
%         commands; |MATPOWER| is installed for this session
%   save_it (integer or string) : indicates whether or not to save the
%       results
%
%       - 0 or [] *(default)* - don't save any results
%       - if ``modify`` is 0
%
%           - ``save_it`` = :samp:`{<string>}` : the relevant :func:`addpath`
%             commands are saved to a file whose name is provided in
%             ``save_it``; execute saved file in any session to make
%             |MATPOWER| available for the session
%           - ``save_it`` = :samp:`{<other true value>}` : the relevant
%             :func:`addpath` commands are saved to a file named
%             ``'startup.m'`` in the current directory; |MATPOWER| is
%             available in any session affected by this ``'startup.m'`` file
%           - *otherwise* : the commands are displayed, but not saved
%       - if ``modify`` is 1
%
%           - ``save_it`` = :samp:`{<any true value>}` : the path will be
%             modified and saved with :func:`savepath`; |MATPOWER| is
%             available in this and all future sessions
%           - *otherwise* : the path will be modified but not saved
%   verbose (boolean) : prints the relevant :func:`addpath` commands if true
%       *(default)*, silent otherwise
%   rm_oldpaths (boolean) : remove existing installation
%
%       - 0 *(default)* - do **not** remove existing |MATPOWER| from path
%       - 1 - remove existing |MATPOWER| paths first
%
% Output:
%     success (boolean) : 1 if all commands succeeded, 0 otherwise
%
% Examples::
%
%   install_matpower                    % interactive mode, prompt for options
%   install_matpower(0);                % print the required addpath() commands
%   install_matpower(0, 1);             % save the commands to startup.m
%   install_matpower(1, 1);             % modify my path and save
%   install_matpower(1, 0, 0);          % modify my path temporarily & silently
%   install_matpower(0, 'matpower8');   % save the commands to matpower8.m
%   install_matpower(0, 0, 1, 1);       % uninstall MATPOWER from path (must
%                                       % call savepath() separately to make
%                                       % permanent)
%
%   See also addpath, savepath.

%   MATPOWER
%   Copyright (c) 2017-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% installation data for each component
min_ver.Octave = '4.0.0';
min_ver.MATLAB = '7.9.0';
min_ver_mp_core.Octave = '6.2.0';
min_ver_mp_core.MATLAB = '9.0.0';
install = struct( ...
    'name', { ...       %% (not currently used)
        'matpower', ...
        'most', ...
        'mp-opt-model', ...
        'mips', ...
        'mptest', ...
        'maxloadlim', ...
        'misc', ...
        'reduction', ...
        'sdp_pf', ...
        'se', ...
        'simulink_matpower', ...
        'smartmarket', ...
        'state_estimator', ...
        'syngrid' }, ...
    'dirs', { ...       %% directories to add to path
        {{'lib'}, {'lib', 't'}, {'data'}}, ...
        {{'most', 'lib'}, {'most', 'lib', 't'}}, ...
        {{'mp-opt-model', 'lib'}, {'mp-opt-model', 'lib', 't'}}, ...
        {{'mips', 'lib'}, {'mips', 'lib', 't'}}, ...
        {{'mptest', 'lib'}, {'mptest', 'lib', 't'}}, ...
        {{'extras', 'maxloadlim'}, {'extras', 'maxloadlim', 'tests'}, ...
            {'extras', 'maxloadlim', 'examples'}}, ...
        {{'extras', 'misc'}}, ...
        {{'extras', 'reduction'}}, ...
        {{'extras', 'sdp_pf'}}, ...
        {{'extras', 'se'}}, ...
        {{'extras', 'simulink_matpower'}, ...
            {'extras', 'simulink_matpower', 'helper_functions'}, ...
            {'extras', 'simulink_matpower', 'helper_functions', 'getter_setter_functions'}}, ...
        {{'extras', 'smartmarket'}}, ...
        {{'extras', 'state_estimator'}}, ...
        {{'extras', 'syngrid', 'lib'}, {'extras', 'syngrid', 'lib', 't'}} }, ...
    'fcns', { ...       %% functions to search for to find old paths to remove
        {'mpver', 'test_matpower', 'case9'}, ...
        {'mostver', 'test_most'}, ...
        {'mpomver', 'test_mp_opt_model'}, ...
        {'mipsver', 'test_mips'}, ...
        {'t_begin', 'test_mptest'}, ...
        {'maxloadlim', 'test_mll_main', 'example_ieee9.m'}, ...
        {'checklimits'}, ...
        {'MPReduction'}, ...
        {'sdp_pf_ver'}, ...
        {'run_se'}, ...
        {'SimulinkMATPOWERbase.slx', 'ac_testbed_setup', 'get_losses_from_matpower'}, ...
        {'runmarket'}, ...
        {'runse'}, ...
        {'syngrid', 'test_syngrid'} } ...
);
ni = length(install);   %% number of components

%% default arguments
interactive = 0;
if nargin < 4
    rm_oldpaths = 0;
    if nargin < 3
        verbose = 1;
        if nargin < 2
            save_it = 0;
            if nargin < 1
                modify = [];
                interactive = 1;
            end
        end
    end
end

%% get path to new MATPOWER installation root
[root, n, e] = fileparts(which('install_matpower'));

%% MATLAB or Octave
if exist('OCTAVE_VERSION', 'builtin') == 5
    sw = 'Octave';
else
    sw = 'MATLAB';
end

%% check for required version of MATLAB or Octave
vstr = '';
switch sw
    case 'Octave'
        v = ver('octave');
    case 'MATLAB'
        v = ver('matlab');
        if length(v) > 1
            warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
            v = v(1);
        end
end
if ~isempty(v) && isfield(v, 'Version') && ~isempty(v.Version)
    vstr = v.Version;
    if vstr2num(vstr) < vstr2num(min_ver.(sw))
        error('\n\n%s\n  MATPOWER requires %s %s or later.\n      You are using %s %s.\n%s\n\n', repmat('!',1,45), sw, min_ver.(sw), sw, vstr, repmat('!',1,45));
    end
else
    warning('\n\n%s\n  Unable to determine your %s version. This indicates\n  a likely problem with your %s installation.\n%s\n', repmat('!',1,60), sw, sw, repmat('!',1,60));
end


%% get installation options interactively, if necessary
div_line = sprintf('\n-------------------------------------------------------------------\n\n');
if interactive
    fprintf(div_line);
    fprintf('MATPOWER Installation Options:\n\n');
    fprintf('   1. Do NOT modify the %s path\n', sw);
    fprintf('          (just generate the required ADDPATH commands)\n');
    fprintf('   2. DO modify the %s path, but only temporarily\n', sw);
    fprintf('          (you will have to do it again next time you run %s)\n', sw);
    fprintf('   3. DO modify the %s path, and SAVE the updated path\n', sw);
    fprintf('          (so you will not have to do it again next time you run %s)\n\n', sw);
    s = 0;
    while ~isempty(s) && (s < 1 || s > 3 || rem(s, 1))
        s = str2num(input('Please enter your selection [1, 2, 3] (default = 1) : ', 's'));
    end
    if isempty(s)
        s = 1;              %% default to option 1
    end
    switch s
        case 1              %% don't modify
            modify = 0;
        case 2              %% modify, don't save
            modify = 1;
            save_it = 0;
        case 3              %% modify and save
            modify = 1;
            save_it = 1;
        otherwise
    end
end

%% find available new components
available = zeros(ni, 1);
for i = 1:ni
    if exist(fullfile(root, install(i).dirs{1}{:}), 'dir')
        available(i) = 1;
    end
    % fprintf('%d %s\n', available(i), fullfile(root, install(i).dirs{:}));
end

%% find paths for currently installed components
oldpaths = {};
for i = 1:ni
    if available(i)     %% only if we have a new replacement
        for k = 1:length(install(i).fcns)
            [p, n, e] = fileparts(which(install(i).fcns{k}));
            if ~isempty(p) && ~ismember(p, oldpaths)
                oldpaths{end+1} = p;
            end
        end
    end
end

if ~isempty(oldpaths)
    if exist('mpver', 'file')
        v = mpver;
    else
        v = '';
    end
    str = sprintf('It appears you already have MATPOWER %s or some of its\n', v);
    str = sprintf('%scomponents installed in the following directories:\n\n', str);
    rm_path_str = '';
    for k = 1:length(oldpaths)
        rm_path_str = sprintf('%s    %s\n', rm_path_str, oldpaths{k});
    end
    if interactive
        fprintf(div_line);
        fprintf('%s%s\n', str, rm_path_str);
        fprintf('Unless you are sure you know what you are doing, we recommend\n');
        fprintf('removing these directories from your %s path before proceeding.\n\n', sw);
        s = '?';
        while ~isempty(s) && (s(1) ~= 'Y' && s(1) ~= 'N')
            s = upper(input('Would you like me to do it for you? [Y or N] (default = N) : ', 's'));
        end
        if isempty(s)
            s = 'N';                %% default to N
        end
        if s(1) == 'Y'
            rm_oldpaths = 1;
        end
    elseif ~rm_oldpaths
        error('install_matpower: %s%s\nPlease remove the old installation first, or re-run\nINSTALL_MATPOWER in interactive mode (no arguments).', str, rm_path_str);
    end

    %% remove old paths
    if rm_oldpaths
        rmpath(oldpaths{:});
        if verbose
            fprintf(div_line);
            fprintf('The following directories were removed from your %s path:\n\n%s\n', sw, rm_path_str);
            if ~modify
                fprintf('You will need to manually use SAVEPATH to make the changes permanent.\n');
                if interactive
                    s = input('Hit any key to continue ...', 's');
                end
            end
        end
    end
end

%% generate paths to add
newpaths = {};
for i = 1:ni
    if available(i)     %% only available components
        for k = 1:length(install(i).dirs)
            p = fullfile(root, install(i).dirs{k}{:});
            if ~isempty(p) && ~ismember(p, newpaths)
                newpaths{end+1} = p;
            end
        end
    end
end

%% build addpath string
cmd = sprintf('addpath( ...\n');
for k = 1:length(newpaths)
    cmd = sprintf('%s    ''%s'', ...\n', cmd, newpaths{k});
end
cmd = sprintf('%s    ''%s'' );\n', cmd, '-end');

%% print addpath
if verbose
    fprintf(div_line);
    if modify
        fprintf('Your %s path will be updated using the following command:\n\n%s', sw, cmd);
        if interactive
            s = input(sprintf('\nHit any key to continue ...'), 's');
        end
    else
        fprintf('Use the following command to add MATPOWER to your %s path:\n\n%s\n', sw, cmd);
    end
end

%% add the new paths
if modify
    addpath(newpaths{:}, '-end');
    if verbose
        fprintf('Your %s path has been updated.\n', sw);
    end
end

%% prompt about saving the command to a file
if interactive && ~modify
    s = '?';
    while ~isempty(s) && (s(1) ~= 'Y' && s(1) ~= 'N')
        s = upper(input('Would you like to save this command to a file? [Y or N] (default = N) : ', 's'));
    end
    if isempty(s)
        s = 'N';                %% default to N
    end
    if s(1) == 'Y'
        %% prompt for file name
        s = '?';
        while ~isempty(s) && s(1) == '?'
            s = input('Please enter the name of file? (default = ''startup.m'') : ', 's');
        end
        if isempty(s)
            s = 'startup.m';    %% default to 'startup.m'
        end
        save_it = s;
    else
        save_it = 0;
    end
end

%% save path
success = 1;
if save_it
    if modify                   %% modify the path directly
        savepath;                   %% save the updated path
        if verbose
            fprintf(div_line);
            fprintf('Your updated %s path has been saved using SAVEPATH.\n', sw);
        end
    else                        %% don't modify the path directly
        %% set up file name
        if ~ischar(save_it)
            save_it = 'startup.m';
        end
        [p, n, e] = fileparts(save_it);
        if isempty(e)
            e = '.m';
        end
        fn = fullfile(p, [n e]);

        %% write file
        if exist(fn, 'file')        %% file exists, print a warning
            if verbose
                fprintf(div_line);
                fprintf('The file ''%s'' was not written. A file with that name already exists.', fn);
            else
                error('install_matpower: file ''%s'' not written, a file with that name already exists', fn);
            end
        else                        %% create the file
            fid = fopen(fn, 'a');
            if fid == -1
                success = 0;
            else
                fprintf(fid, '%%%s\n\n', upper(n));
                fprintf(fid, '%%%% add MATPOWER paths\n');
                fprintf(fid, '%s', cmd);
                fclose(fid);
            end
            if verbose
                fprintf(div_line);
                fprintf('The file ''%s'' containing the commands to\nadd MATPOWER to your %s path has been created.\n', fn, sw);
            end
        end
    end
end

if verbose
    fprintf(div_line);
    if modify
        fprintf('Now that you have added the required directories to your %s path\n', sw);
        fprintf('MATPOWER is installed and ready to use.\n\n');
    else
        fprintf('Once you have added the required directories to your %s path\n', sw);
        fprintf('MATPOWER will be installed and ready to use.\n\n');
    end
    fprintf('You may want to begin by typing: mpver\n');
    fprintf('to see the list of installed MATPOWER related software versions.\n\n');
    fprintf('Or to run the MATPOWER test suite to ensure everything is\n');
    fprintf('working correctly, type: test_matpower\n\n');

    if interactive && modify
        s = '?';
        while ~isempty(s) && (s(1) ~= 'Y' && s(1) ~= 'N')
            s = upper(input('Would you like to run the MATPOWER tests now? [Y or N] (default = N) : ', 's'));
        end
        if isempty(s)
            s = 'N';                %% default to N
        end
        if s(1) == 'Y'
            test_matpower;
        end
    end
end

if nargout
    succ = success;
end

function num = vstr2num(vstr)
% Converts version string to numerical value suitable for < or > comparisons
% E.g. '3.11.4' -->  3.011004
pat = '\.?(\d+)';
[s,e,tE,m,t] = regexp(vstr, pat);
b = 1;
num = 0;
for k = 1:length(t)
    num = num + b * str2num(t{k}{1});
    b = b / 1000;
end
