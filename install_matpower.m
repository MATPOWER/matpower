function succ = install_matpower(modify, save_it, verbose)
%INSTALL_MATPOWER  Assist the user in setting their path.
%   INSTALL_MATPOWER
%   INSTALL_MATPOWER(MODIFY)
%   INSTALL_MATPOWER(MODIFY, SAVE_IT)
%   INSTALL_MATPOWER(MODIFY, SAVE_IT, VERBOSE)
%   SUCCESS = INSTALL_MATPOWER(...)
%
%   Assists the user in setting up the proper MATLAB/Octave path to
%   be able to use MATPOWER and run its tests. With no input arguments
%   it prompts interactively to determine how to handle the paths.
%
%   Inputs (all are optional):
%       MODIFY : select how to set path
%           0 (default) - generate relevant ADDPATH() commands, but
%               don't execute them
%           1 - modify the path directly executing the relevant
%               ADDPATH() commands
%       SAVE_IT : indicates whether or not to save the results
%           0 or [] (default) - don't save any results
%           if MODIFY is 0
%               SAVE_IT = <string> : the relevant ADDPATH() commands
%                   are saved to a file whose name is provided in SAVE_IT
%               SAVE_IT = <other true value> : the relevant ADDPATH() commands
%                   are saved to a file named 'startup.m' in the current
%                   directory
%               otherwise : the commands are displayed, but not saved
%           if MODIFY is 1
%               SAVE_IT = <any true value> : the path will be modified
%                   and saved with SAVEPATH
%               otherwise : the path will be modified but not saved
%       VERBOSE : prints the relevant ADDPATH commands if true (default),
%               silent otherwise
%
%   Outputs (all are optional):
%       SUCCESS : 1 if all commands succeeded, 0 otherwise
%
%   Examples:
%       install_matpower;           %% print the required ADDPATH() commands
%       install_matpower(0, 1);     %% save the commands to startup.m
%       install_matpower(1, 1);     %% modify my path and save
%       install_matpower(1, 0, 0);  %% modify my path temporarily and silently
%       install_matpower(0, 'matpower6');   %% save the commands to matpower6.m
%
%   See also ADDPATH, SAVEPATH.

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% installation data for each component
min_ver.Octave = '4.0.0';
min_ver.MATLAB = '7.3.0';
install = struct( ...
    'name', { ...
        'matpower', ...
        'mips', ...
        'most', ...
        'mptest', ...
        'maxloadlim', ...
        'misc', ...
        'reduction', ...
        'sdp_pf', ...
        'se', ...
        'smartmarket', ...
        'state_estimator' }, ...
    'dirs', { ...
        {{'lib'}, {'lib', 't'}, {'data'}}, ...
        {{'mips', 'lib'}, {'mips', 'lib', 't'}}, ...
        {{'most', 'lib'}, {'most', 'lib', 't'}}, ...
        {{'mptest', 'lib'}, {'mptest', 'lib', 't'}}, ...
        {{'extras', 'maxloadlim'}, {'extras', 'maxloadlim', 'examples'}}, ...
        {{'extras', 'misc'}}, ...
        {{'extras', 'reduction'}}, ...
        {{'extras', 'sdp_pf'}}, ...
        {{'extras', 'se'}}, ...
        {{'extras', 'smartmarket'}}, ...
        {{'extras', 'state_estimator'}} }, ...
    'fcns', { ...
        {'mpver', 'test_matpower', 'case9'}, ...
        {'mipsver', 'test_mips'}, ...
        {'mostver', 'test_most'}, ...
        {'t_begin', 'test_mptest'}, ...
        {'maxloadlim', 'example_ieee9.m'}, ...
        {'checklimits'}, ...
        {'MPReduction'}, ...
        {'sdp_pf_ver'}, ...
        {'run_se'}, ...
        {'runmarket'}, ...
        {'runse'} } ...
);
ni = length(install);   %% number of components

%% default arguments
interactive = 0;
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

rm_oldpaths = 0;
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
    else    
        error('install_matpower: %s%s\nPlease remove the old installation first, or re-run\nINSTALL_MATPOWER in interactive mode (no arguments).', str, rm_path_str);
    end
end

%% remove old paths
if rm_oldpaths
    rmpath(oldpaths{:});
    if verbose
        fprintf(div_line);
        fprintf('The following directories were removed from your %s path:\n\n%s\n', sw, rm_path_str);
        if ~modify
            fprintf('You will need to manually use SAVEPATH to make the changes permanent.\n');
            s = input('Hit any key to continue ...', 's');
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
