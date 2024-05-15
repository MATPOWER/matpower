function [success, err_list] = generate_autodoc_stubs(in, dest)
% generate_autodoc_stubs - Generate source .rst stubs for Reference Manual.
% ::
%
%   generate_autodoc_stubs(in, dest)
%
%   % not yet implemented
%   [success, err_list] = generate_autodoc_stubs(in, dest)
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
%       generate_autodoc_stubs(in, ...
%           '~/matpower/docs/sphinx/source/ref-manual/');

%   MATPOWER
%   Copyright (c) 2023-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% check for destination directory
if ~exist(dest, 'dir')
    error('generate_autodoc_stubs: directory ''%s'' does not exist.', dest);
end

stub_types = fieldnames(in);
for k = 1:length(stub_types)
    stub_type = in.(stub_types{k});
    fprintf('\nCreating .rst stub files for %s:\n', stub_types{k});

    %% create destination directory
    dest_dir = fullfile(dest, stub_type.destdir);
    if ~exist(dest_dir, 'dir')
        fprintf('Creating directory ''%s''\n', dest_dir);
        success = mkdir(dest, stub_type.destdir);
        if ~success
            error('generate_autodoc_stubs: unable to create directory ''%s'' in ''%s''.', ...
                stub_type.destdir, dest);
        end
    end

    %% create stub files
    for m = 1:length(stub_type.list)
        mod = stub_type.list(m).mod;
        rel_paths = strsplit(mod, '.');
        rel_paths = rel_paths(2:end);
        pkg = strrep(mod, '+', '');     %% strip any plus signs from pkg
        pkgs = strsplit(pkg, '.');
        pkgs = pkgs(2:end);             %% eliminate module directory name
        pkg_dir = fullfile(dest_dir, pkgs{:});
        if ~exist(pkg_dir, 'dir')
            fprintf('Creating directory ''%s''\n', pkg_dir);
            success = mkdir(pkg_dir);
            if ~success
                error('generate_autodoc_stubs: unable to create directory ''%s''.', ...
                    pkg_dir);
            end
        end
        names = stub_type.list(m).names;
        for f = 1:length(names)
            %% open file
            if isempty(pkgs)
                pkg_str = '';
            else
                pkg_str = [strjoin(pkgs, '/') '/'];
            end
            if names{f}(1) == '@'   %% class folder
                is_class_folder = 1;
                names{f} = names{f}(2:end); %% strip the @
            else
                is_class_folder = 0;
            end

            fprintf('    %s%s.rst\n', pkg_str, names{f});
            fname = fullfile(pkg_dir, [names{f} '.rst']);
            fullname = strjoin({pkgs{:}, names{f}}, '.');
            if is_class_folder
                relpath = strjoin({rel_paths{:}, ['@' names{f}], names{f}}, '/');
            else
                relpath = strjoin({rel_paths{:}, names{f}}, '/');
            end
            [fd, msg] = fopen(fname, 'wt');
%             fprintf(fd, '.. toctree::\n');
%             fprintf(fd, '   :maxdepth: 2\n');
%             fprintf(fd, '\n');
            fprintf(fd, '.. currentmodule:: %s\n', mod);
            fprintf(fd, '\n');
            if isempty(stub_type.list(m).src_path)
                stub_path = relpath;
            else
                stub_path = sprintf('%s/%s', stub_type.list(m).src_path, relpath);
            end
            fprintf(fd, ':raw-html:`<div style="float: right"><a href="%s/%s.m" target=_blank><svg height="32" aria-hidden="true" viewBox="0 0 16 16" version="1.1" width="32" data-view-component="true" class="octicon octicon-mark-github v-align-middle color-fg-default"><path d="M8 0c4.42 0 8 3.58 8 8a8.013 8.013 0 0 1-5.45 7.59c-.4.08-.55-.17-.55-.38 0-.27.01-1.13.01-2.2 0-.75-.25-1.23-.54-1.48 1.78-.2 3.65-.88 3.65-3.95 0-.88-.31-1.59-.82-2.15.08-.2.36-1.02-.08-2.12 0 0-.67-.22-2.2.82-.64-.18-1.32-.27-2-.27-.68 0-1.36.09-2 .27-1.53-1.03-2.2-.82-2.2-.82-.44 1.1-.16 1.92-.08 2.12-.51.56-.82 1.28-.82 2.15 0 3.06 1.86 3.75 3.64 3.95-.23.2-.44.55-.51 1.07-.46.21-1.61.55-2.33-.66-.15-.24-.6-.83-1.23-.82-.67.01-.27.38.01.53.34.19.73.9.82 1.13.16.45.68 1.31 2.69.94 0 .67.01 1.3.01 1.49 0 .21-.15.45-.55.38A7.995 7.995 0 0 1 0 8c0-4.42 3.58-8 8-8Z"></path></svg></a></div>`\n', ...
                stub_type.gh_base_url, stub_path);
            fprintf(fd, '\n');
            fprintf(fd, '%s\n', fullname);
            fprintf(fd, '%s\n', repmat('-', 1, length(fullname)));
            fprintf(fd, '\n');
            if is_class_folder
                fprintf(fd, '.. auto%s:: %s.@%s.%s\n', stub_types{k}, mod, names{f}, names{f});
            else
                fprintf(fd, '.. auto%s:: %s\n', stub_types{k}, names{f});
            end
            switch stub_types{k}
            case 'class'
                fprintf(fd, '    :show-inheritance:\n');
                fprintf(fd, '    :members:\n');
            case 'function'
            end
            fclose(fd);
        end
    end
end
success = 1;
err_list = {};
