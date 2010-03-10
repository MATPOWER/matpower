function rv = mpver(varargin)
%MPVER  Prints or returns MATPOWER version info for current installation.
%   V = MPVER returns the current MATPOWER version number.
%   V = MPVER('all') returns a struct with the fields Name, Version,
%   Release and Date (all strings). Calling MPVER without assigning the
%   return value prints the version and release date of the current
%   installation of MATPOWER, MATLAB, the Optimization Toolbox, MIPS
%   and any optional MATPOWER packages such as BPMPD_MEX, MINOPF,
%   PDIPMOPF, SCPDIPMOPF and TRAMLOPF.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005-2010 by Power System Engineering Research Center (PSERC)
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
                'Version',  '4.0-dev', ...
                'Release',  '', ...
                'Date',     '03-Mar-2010' );
if nargout > 0
    if nargin > 0
        rv = v{1};
    else
        rv = v{1}.Version;
    end
else
    v{2} = ver('matlab');
    v{3} = ver('optim');
    for n = 1:3
        if n == 3 && isempty(v{3})
            fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
            continue;
        end
        fprintf('\n%-22s Version %-9s  %11s', v{n}.Name, v{n}.Version, v{n}.Date);
        if ~isempty(v{n}.Release)
            fprintf('   Release: %-10s', v{n}.Release);
        end
    end
    fprintf('\n');
    mipsver;
    if have_fcn('bpmpd')
        if exist('bpver', 'file') == 2
            bpver;
        else
            fprintf('BPMPD_MEX              Version 2.21 or earlier\n');
        end
    else
        fprintf('BPMPD_MEX              not installed\n');
    end
    
    if have_fcn('minopf')
        if exist('minopfver', 'file') == 2
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
