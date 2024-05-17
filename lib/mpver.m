function rv = mpver(varargin)
% mpver - Prints or returns installed |MATPOWER| version info.
% ::
%
%   mpver
%   v = mpver
%   v = mpver('all')
%
% When called with an output argument and no input argument, mpver
% returns the current |MATPOWER| version numbers. With an input argument
% (e.g. ``'all'``) it returns  a struct with the fields ``Name``,
% ``Version``, ``Release``, and ``Date`` *(all char arrays)*. Calling
% mpver without assigning the return value prints the version and
% release date of the current installation of |MATPOWER|, |MATLAB| (or
% |OCTAVE|), the Optimization Toolbox, MP-Test, MIPS, MP-Opt-Model, MOST,
% and any optional |MATPOWER| packages.

%   MATPOWER
%   Copyright (c) 2005-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

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
                'Version',  '8.0', ...
                'Release',  '', ...
                'Date',     '17-May-2024' );
if nargout > 0
    if nargin > 0
        rv = v{1};
    else
        rv = v{1}.Version;
    end
else
    if have_feature('octave')
        % v{2} = ver('octave');
        %% workaround for https://savannah.gnu.org/bugs/index.php?59125
        vv = ver;
        for k = 1:length(vv)
            if strcmp(vv(k).Name, 'Octave')
                v{2} = vv(k);
                break;
            end
        end
        %% convert to consistent format
        if ~isempty(v{2}.Date)
            v{2}.Date = datestr(v{2}.Date, 'dd-mmm-yyyy');
        end
    else
        v{2} = ver('matlab');
        if length(v{2}) > 1
            warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''matlab'' on your path. Check each element of the output of ver(''matlab'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
            v{2} = v{2}(1);
        end
    end
    v{3} = ver('optim');
    if length(v{3}) > 1
        warning('The built-in VER command is behaving strangely, probably as a result of installing a 3rd party toolbox in a directory named ''optim'' on your path. Check each element of the output of ver(''optim'') to find the offending toolbox, then move the toolbox to a more appropriately named directory.');
        v{3} = v{3}(1);
    end
    for n = 1:3
        if n == 3
            if isempty(v{3}) || isempty(v{3}.Version)
                if have_feature('matlab')
                    fprintf('\n%-22s -- not installed --', 'Optimization Toolbox');
                else    %% Octave
                    fprintf('\n%-22s -- not installed --', 'optim');
                end
                continue;
            elseif have_feature('matlab') && ~license('test', 'optimization_toolbox')
                fprintf('\n%-22s -- no license --', 'Optimization Toolbox');
                continue;
            end
        end
        fprintf('\n%-22s Version %-9s', v{n}.Name, v{n}.Version);
        if ~isempty(v{n}.Date)
            fprintf('  %11s', v{n}.Date);
            if ~isempty(v{n}.Release)
                fprintf('   Release: %-10s', v{n}.Release);
            end
        end
    end
    fprintf('\n');
    mptestver;
    mipsver;
    mpomver;
    if have_feature('most')
        mostver;
    else
        fprintf('%-22s -- not installed --\n', 'MOST');
    end
    if have_feature('e4st')
        e4st_ver;
    end
    if have_feature('sdp_pf')
        sdp_pf_ver;
    else
        fprintf('%-22s -- not installed --\n', 'SDP_PF');
    end
    if have_feature('syngrid')
        sgver;
    else
        fprintf('%-22s -- not installed --\n', 'SynGrid');
    end
    if have_feature('knitro')
        s = have_feature('knitro', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'Artelys Knitro', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'Artelys Knitro');
    end
    if have_feature('bpmpd')
        if exist('bpver', 'file') == 2
            bpver;
        else
            fprintf('%-22s Version 2.21 or earlier\n', 'BPMPD_MEX');
        end
    else
        fprintf('%-22s -- not installed --\n', 'BPMPD_MEX');
    end
    if have_feature('clp')
        s = have_feature('clp', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'CLP', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CLP');
    end
    if have_feature('cplex')
        s = have_feature('cplex', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'CPLEX', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'CPLEX');
    end
    if have_feature('glpk')
        s = have_feature('glpk', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'GLPK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'GLPK');
    end
    gurobiver;
    if have_feature('ipopt')
        s = have_feature('ipopt', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'IPOPT', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'IPOPT');
    end
    if have_feature('minopf')
        if exist('minopfver', 'file') == 2
            minopfver;
        else
            fprintf('%-22s Version 3.0b2 or earlier\n', 'MINOPF');
        end
    else
        fprintf('%-22s -- not installed --\n', 'MINOPF');
    end
    if have_feature('mosek')
        s = have_feature('mosek', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'MOSEK', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'MOSEK');
    end
    osqpver;
    if have_feature('pardiso')
        s = have_feature('pardiso', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'PARDISO', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'PARDISO');
    end
    if have_feature('pdipmopf')
        pdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'PDIPMOPF');
    end
    if have_feature('scpdipmopf')
        scpdipmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'SCPDIPMOPF');
    end
    if have_feature('sdpt3')
        s = have_feature('sdpt3', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SDPT3', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SDPT3');
    end
    if have_feature('sedumi')
        s = have_feature('sedumi', 'all');
        if isempty(s.vstr)
            vn = '<unknown>';
        else
            vn = s.vstr;
        end
        fprintf('%-22s Version %-10s %-11s\n', 'SeDuMi', vn, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'SeDuMi');
    end
    if have_feature('tralmopf')
        tralmopfver;
    else
        fprintf('%-22s -- not installed --\n', 'TRALMOPF');
    end
    if have_feature('yalmip')
        s = have_feature('yalmip', 'all');
        fprintf('%-22s Version %-10s %-11s\n', 'YALMIP', s.vstr, s.date);
    else
        fprintf('%-22s -- not installed --\n', 'YALMIP');
    end

    fprintf('%-22s %s\n\n', 'Architecture:', computer);
    
    fprintf('  MATPOWER %s is distributed under the 3-clause BSD License.\n', v{1}.Version);
    fprintf('  Please see the LICENSE file for details.\n\n');
end
