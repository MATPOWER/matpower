function [success, err_list] = generate_source_symlinks(in, dest, src)
% generate_source_symlinks - Generate source symlinks for Reference Manual.
% ::
%
%   generate_source_symlinks(in, dest, src)
%
%   % not yet implemented
%   [success, err_list] = generate_source_symlinks(in, dest, src)
%
% Inputs:
%   in (struct) : input specification for source files of interest, with fields:
%       ::
%
%           in
%               .(stub_type)    e.g. 'class', 'function'
%               .destdir        prepend dest to get dir for .rst stub files,
%                               e.g. 'classes', 'functions'
%               .gh_base_url    base URL for GitHub source
%               .list           struct array
%                   .mod        module names for Sphinx, e.g. 'matpower',
%                               'matpower.+pkg'
%                   .src_path   relative path from base of source to this
%                               set of source files (e.g. 'lib', 'lib/t')
%                   .names      cell array of function/class names
%                               (w/o .m extension)
%
%   dest (char array) : path the base directory for symlinks
%   src (char array) : relative path of the base directory of the original
%       source from ``dest``
%
% Outputs:
%   success (integer) : *(not yet implemented)* - 1 if successful, 0 otherwise
%   err_list (cell array) : *(not yet implemented)* - all errors are currently
%       fatal
%
% Example::
%
%       in = struct(...
%           'class', struct(...
%               'destdir', 'classes', ...
%               'gh_base_url', 'https://github.com/MATPOWER/matpower/blob/master', ...
%               'list', struct(...
%                   'mod', {'matpower', 'matpower.+mp'}, ...
%                   'src_path', {'lib', 'lib/t'}, ...
%                   'names', {{'class1', 'class2'}, {'class3', 'class4'}} ... - 
%               ) ...
%           ), ...
%           'function', struct(...
%               'destdir', 'functions', ...
%               'gh_base_url', 'https://github.com/MATPOWER/matpower/blob/master', ...
%               'list', struct(...
%                   'mod', {'matpower', 'matpower.+mp'}, ...
%                   'src_path', {'lib', 'lib'}, ...
%                   'names', {{'func1', 'func2'}, {'func3', 'func4'}} ... - 
%               ) ...
%           ) ...
%       );
%       generate_source_symlinks(in, ...
%           '~/matpower/docs/sphinx/source/matlab-source/', ...
%           '../../../../../');

%   MATPOWER
%   Copyright (c) 2023-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% check for destination directory
if ~exist(dest, 'dir')
    error('generate_source_symlinks: directory ''%s'' does not exist.', dest);
end

cwd = pwd;
stub_types = fieldnames(in);
for k = 1:length(stub_types)
    stub_type = in.(stub_types{k});
    fprintf('\nCreating symlinks for %s:\n', stub_types{k});

    %% create symlinkss
    for m = 1:length(stub_type.list)
        mod = stub_type.list(m).mod;
        rel_paths = strsplit(mod, '.');

        %% create destination directory if necessary
        dest_dir = fullfile(dest, rel_paths{:});
        if ~exist(dest_dir, 'dir')
            fprintf('Creating directory ''%s''\n', dest_dir);
            success = mkdir(dest_dir);
            if ~success
                error('generate_source_symlinks: unable to create directory ''%s''.', ...
                    destdir);
            end
        end

        %% change directory to dest_dir
        cd(dest_dir);

        rel_paths = rel_paths(2:end);
        pkg_dir = strjoin(rel_paths, '/');
        if ~isempty(pkg_dir)
            pkg_dir = ['/' pkg_dir];
        end
        rel_path_str = repmat('../', 1, length(rel_paths));
        tgt_dir = sprintf('%s%s%s%s', src, rel_path_str, stub_type.list(m).src_path, pkg_dir);
        if tgt_dir(end) == '/'
            tgt_dir(end) = [];  % delete trailing '/' (e.g. empty src_dir)
        end
        names = stub_type.list(m).names;
        for f = 1:length(names)
            %% create symlink
            if names{f}(1) == '@'
                sl_src = sprintf('%s', names{f});
                sl_tgt = sprintf('%s/%s', tgt_dir, names{f});
            else
                sl_src = sprintf('%s.m', names{f});
                sl_tgt = sprintf('%s/%s.m', tgt_dir, names{f});
            end
            if exist(sl_tgt)
                fprintf('    %s/%s\n', stub_type.list(m).src_path, sl_src);
                if exist(['./' sl_src])
                    delete(['./' sl_src]);
                end
                cmd = sprintf('ln -s %s %s', sl_tgt, sl_src);
                system(cmd);
            else
                fprintf('    ----> target for %s does NOT exist\n', sl_src);
            end
            
        end
    end
end
cd(cwd);
success = 1;
err_list = {};
