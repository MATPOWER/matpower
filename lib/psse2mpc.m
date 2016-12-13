function [mpc, warnings] = psse2mpc(rawfile_name, mpc_name, verbose, rev)
%PSSE2MPC  Converts a PSS/E RAW data file into a MATPOWER case struct.
%   MPC = PSSE2MPC(RAWFILE_NAME)
%   MPC = PSSE2MPC(RAWFILE_NAME, VERBOSE)
%   MPC = PSSE2MPC(RAWFILE_NAME, VERBOSE, REV)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME, VERBOSE)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME, VERBOSE, REV)
%   [MPC, WARNINGS] = PSSE2MPC(RAWFILE_NAME, ...)
%
%   Converts a PSS/E RAW data file into a MATPOWER case struct.
%
%   Input:
%       RAWFILE_NAME :  name of the PSS/E RAW file to be converted
%                       (opened directly with FILEREAD)
%       MPC_NAME     :  (optional) file name to use to save the resulting
%                        MATPOWER case
%       VERBOSE      :  1 (default) to display progress info, 0 otherwise
%       REV          :  (optional) assume the input file is of this
%                       PSS/E revision number, attempts to determine
%                       REV from the file by default
%
%   Output(s):
%       MPC      : resulting MATPOWER case struct
%       WARNINGS : (optional) cell array of strings containing warning
%                  messages (included by default in comments of MPC_NAME).
%
% NOTE: The data sections to be read in the PSS/E raw file includes:
%       identification data; bus data; branch data; fixed shunt data;
%       generator data; transformer data; switched shunt data; area data
%       and hvdc line data. Other data sections are currently ignored.

%   MATPOWER
%   Copyright (c) 2014-2016, Power Systems Engineering Research Center (PSERC)
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mpraw2mp.m, written by: Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% handle input args
if nargin < 2
    rev = 0;
    verbose = 1;
    mpc_name = '';
elseif ischar(mpc_name)     %% save the file
    if nargin < 4
        rev = 0;
        if nargin < 3
            verbose = 1;
        end
    elseif isempty(rev)
        rev = 0;
    end
else                        %% don't save the file
    if nargin < 3
        rev = 0;
    else
        rev = verbose;
    end
    verbose = mpc_name;
    mpc_name = '';
end

%% read data from PSS/E RAW file
[records, sections] = psse_read(rawfile_name, verbose);

%% parse data
[data, warnings] = psse_parse(records, sections, verbose, rev);

%% convert to MATPOWER case file
[mpc, warnings] = psse_convert(warnings, data, verbose);

%% (optionally) save MATPOWER case file
if ~isempty(mpc_name)
    if ~rev
        rev = data.id.REV;
    end
    comments = {''};
    for k = 0:2
        str = data.id.(sprintf('comment%d', k));
        if ~isempty(str)
            comments{end+1} = sprintf('   %s', str);
        end
    end
    comments{end+1} = '';
    comments{end+1} = sprintf('   Converted by MATPOWER %s using PSSE2MPC on %s', mpver, date);
    comments{end+1} = sprintf('   from ''%s'' using PSS/E rev %d format.', rawfile_name, rev);

    %% warnings
    comments{end+1} = '';
    comments{end+1} = '   WARNINGS:';
    for k = 1:length(warnings)
        comments{end+1} = sprintf('       %s', warnings{k});
    end
    comments{end+1} = '';
    comments{end+1} = sprintf('   See CASEFORMAT for details on the MATPOWER case file format.');

    if verbose
        spacers = repmat('.', 1, 45-length(mpc_name));
        fprintf('Saving to MATPOWER case ''%s'' %s', mpc_name, spacers);
    end
    savecase(mpc_name, comments, mpc);
    if verbose
        fprintf(' done.\n');
    end
end
