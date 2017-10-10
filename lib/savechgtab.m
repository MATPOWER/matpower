function fname_out = savechgtab(fname, chgtab, warnings)
%SAVECHGTAB  Save a change table to a file.
%   SAVECHGTAB(FNAME, CHGTAB)
%   SAVECHGTAB(FNAME, CHGTAB, WARNINGS)
%   FNAME = SAVECHGTAB(FNAME, ...)
%
%   Writes a CHGTAB, suitable for use with APPLY_CHANGES to a file
%   specified by FNAME. If FNAME ends with '.mat' it saves CHGTAB
%   and WARNINGS to a MAT-file as the variables 'chgtab' and 'warnings',
%   respectively. Otherwise, it saves an M-file function that returns
%   the CHGTAB, with the optional WARNINGS in comments.
%
%   Optionally returns the filename, with extension added if necessary.
%
%   Input:
%       FNAME :  name of the file to be saved
%       CHGTAB :  change table suitable for use with APPLY_CHANGES
%       WARNINGS : optional cell array of warning messages (to be
%                   included in comments), such as those returned by
%                   PSSECON2CHGTAB
%
%   Output(s):
%       FNAME :  name of the file, with extention added if necessary

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% handle input args
if nargin < 3
    warnings = {};
end

%% define named indices into data matrices
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% constant names
bus_cols = {'BUS_I', 'BUS_TYPE', 'PD', 'QD', 'GS', 'BS', 'BUS_AREA', ...
    'VM', 'VA', 'BASE_KV', 'ZONE', 'VMAX', 'VMIN'};
br_cols = {'F_BUS', 'T_BUS', 'BR_R', 'BR_X', 'BR_B', 'RATE_A', ...
    'RATE_B', 'RATE_C', 'TAP', 'SHIFT', 'BR_STATUS', 'ANGMIN', 'ANGMAX'};
gen_cols = {'GEN_BUS', 'PG', 'QG', 'QMAX', 'QMIN', 'VG', 'MBASE', ...
    'GEN_STATUS', 'PMAX', 'PMIN', 'PC1', 'PC2', 'QC1MIN', 'QC1MAX', ...
    'QC2MIN', 'QC2MAX', 'RAMP_AGC', 'RAMP_10', 'RAMP_30', 'RAMP_Q', 'APF'};
gc_cols = {'CT_MODCOST_X', 'CT_MODCOST_F', '0', ...  %% use val+3
    'MODEL', 'STARTUP', 'SHUTDOWN', 'NCOST', 'COST'};
load_cols = {'CT_LOAD_ALL_PQ', 'CT_LOAD_FIX_PQ', 'CT_LOAD_DIS_PQ', ...
    'CT_LOAD_ALL_P', 'CT_LOAD_FIX_P', 'CT_LOAD_DIS_P'};
chgtypes = {'CT_REP', 'CT_REL', 'CT_ADD'};

%% verify valid filename
[pathstr, fcn_name, extension] = fileparts(fname);
if isempty(extension)
    extension = '.m';
end
if regexp(fcn_name, '\W')
    old_fcn_name = fcn_name;
    fcn_name = regexprep(fcn_name, '\W', '_');
    fprintf('WARNING: ''%s'' is not a valid function name, changed to ''%s''\n', old_fcn_name, fcn_name);
end
fname = fullfile(pathstr, [fcn_name extension]);

%% open and write the file
if strcmp(upper(extension), '.MAT')     %% MAT-file
    save(fname, 'chgtab', 'warnings');
else                                %% M-file
    %% open file
    [fd, msg] = fopen(fname, 'wt');     %% print it to an M-file
    if fd == -1
        error(['savechgtab: ', msg]);
    end
    
    %% function header, etc.
    fprintf(fd, 'function chgtab = %s\n', fcn_name);
    fprintf(fd, '%%%s  Returns a change table suitable for use with APPLY_CHANGES\n', upper(fcn_name));

    if ~isempty(warnings)
        fprintf(fd, '%%\n');
        fprintf(fd, '%% Warnings:\n');
        for k = 1:length(warnings)
            fprintf(fd, '%%     %s\n', warnings{k});
        end
    end

    fprintf(fd, '\ndefine_constants;\n');

    fprintf(fd, '\n%%%% Change Table\n');
    fprintf(fd, '%%\tlabel\tprob\ttable\trow\tcol\tchgtype\tnewval\n');
    fprintf(fd, 'chgtab = [\n');
    for c = 1:size(chgtab, 1)
        fprintf(fd, '\t%d\t%g\t', ...       %% CT_LABEL, CT_PROB
            chgtab(c, CT_LABEL), chgtab(c, CT_PROB)');
        switch chgtab(c, CT_TABLE)          %% CT_TABLE
            case CT_TBUS
                fprintf(fd, '%s\t', 'CT_TBUS');
                col = bus_cols{chgtab(c, CT_COL)};
            case CT_TGEN
                fprintf(fd, '%s\t', 'CT_TGEN');
                col = gen_cols{chgtab(c, CT_COL)};
            case CT_TGENCOST
                fprintf(fd, '%s\t', 'CT_TGENCOST');
                if chgtab(c, CT_COL)+3 > 8
                    col = sprintf('%s+%d', gc_cols{5+3}, chgtab(c, CT_COL)-5);
                else
                    col = gc_cols{chgtab(c, CT_COL)+3};
                end
            case CT_TBRCH
                fprintf(fd, '%s\t', 'CT_TBRCH');
                col = br_cols{chgtab(c, CT_COL)};
            case CT_TLOAD
                fprintf(fd, '%s\t', 'CT_TLOAD');
                col = load_cols{abs(chgtab(c, CT_COL))};
                if chgtab(c, CT_COL) < 0
                    col = ['-' col];
                end
            case CT_TAREABUS
                fprintf(fd, '%s\t', 'CT_TAREABUS');
                col = bus_cols{chgtab(c, CT_COL)};
            case CT_TAREAGEN
                fprintf(fd, '%s\t', 'CT_TAREAGEN');
                col = gen_cols{chgtab(c, CT_COL)};
            case CT_TAREAGENCOST
                fprintf(fd, '%s\t', 'CT_TAREAGENCOST');
                if chgtab(c, CT_COL)+3 > 8
                    col = sprintf('%s+%d', gc_cols{5+3}, chgtab(c, CT_COL)-5);
                else
                    col = gc_cols{chgtab(c, CT_COL)+3};
                end
            case CT_TAREABRCH
                fprintf(fd, '%s\t', 'CT_TAREABRCH');
                col = br_cols{chgtab(c, CT_COL)};
            case CT_TAREALOAD
                fprintf(fd, '%s\t', 'CT_TAREALOAD');
                col = load_cols{abs(chgtab(c, CT_COL))};
                if chgtab(c, CT_COL) < 0
                    col = ['-' col];
                end
        end
        fprintf(fd, '%g\t%s\t%s\t%g;\n', ...%% CT_ROW, %% CT_COL, CT_CHGTYPE, CT_NEWVAL
            chgtab(c, CT_ROW), col, chgtypes{chgtab(c, CT_CHGTYPE)}, chgtab(c, CT_NEWVAL));
    end
    fprintf(fd, '];\n');

    %% close file
    if fd ~= 1
        fclose(fd);
    end
end

if nargout > 0
    fname_out = fname;
end
