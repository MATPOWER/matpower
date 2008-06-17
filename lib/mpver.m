function mpv = mpver
%MPVER  Prints or returns MATPOWER version information
%   v = mpver returns the current MATPOWER version number. Calling mpver
%   without assigning the return value prints the version and release date of
%   the current installation of MATPOWER, MATLAB, the Optimization Toolbox,
%   and any optional MATPOWER packages such as BPMPD_MEX and MINOPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% the following only works if MATPOWER is explicitly in the path,
%% but not if it is only in the current working directory
% fs = filesep;
% p = fileparts(which('runpf'));
% if ~strcmp(p(1),fs)
%   [t, p] = strtok(p, filesep);
% end
% p = p(2:end);
% v{1} = ver(p);

v{1} = struct(  'Name',     'MATPOWER', ... 
                'Version',  '3.3-dev', ...
                'Release',  '', ...
                'Date',     '17-Jun-2008' );
v{2} = ver('matlab');
v{3} = ver('optim');
if nargout > 0
    mpv = v{1}.Version;
else
    for n = 1:3
        if n == 3 & isempty(v{3})
            fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
            continue;
        end
        fprintf('\n%-22s Version %-9s  %11s', v{n}.Name, v{n}.Version, v{n}.Date);
        if ~isempty(v{n}.Release)
            fprintf('   Release: %-10s', v{n}.Release);
        end
    end
    fprintf('\n');
    
    if have_fcn('bpmpd')
        if exist('bpver')
            bpver
        else
            fprintf('BPMPD_MEX              Version 2.21 or earlier\n');
        end
    else
        fprintf('BPMPD_MEX              not installed\n');
    end
    
    if have_fcn('minopf')
        if exist('minopfver')
            minopfver;
        else
            fprintf('MINOPF                 Version 3.0b2 or earlier\n');
        end
    else
        fprintf('MINOPF                 not installed\n');
    end
    if have_fcn('pdipmopf')
        pdipmopfver;
    else
        fprintf('PDIPMOPF               not installed\n');
    end
    if have_fcn('scpdipmopf')
        scpdipmopfver;
    else
        fprintf('SCPDIPMOPF             not installed\n');
    end
    if have_fcn('tralmopf')
        tralmopfver;
    else
        fprintf('TRALMOPF               not installed\n');
    end

    fprintf('Architecture:          %s\n\n', computer);
end

return;
