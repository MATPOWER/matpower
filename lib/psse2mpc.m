function mpc = psse2mpc(rawfile_name, mpc_name, verbose, rev)
%PSSE2MPC  Converts a PSS/E RAW data file into a MATPOWER case struct.
%   MPC = PSSE2MPC(RAWFILE_NAME)
%   MPC = PSSE2MPC(RAWFILE_NAME, VERBOSE)
%   MPC = PSSE2MPC(RAWFILE_NAME, VERBOSE, REV)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME, VERBOSE)
%   MPC = PSSE2MPC(RAWFILE_NAME, MPC_NAME, VERBOSE, REV)
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
%   Output:
%       MPC : resulting MATPOWER case struct
%
% NOTE: the data sections to be read in the PSS/E raw file includes:
%       identification data; bus data; branch data; fixed shunt data;
%       generator data; transformer data; switched shunt data; area data
%       and hvdc line data
%       other data sections are ignored

%   MATPOWER
%   $Id$
%   by Yujia Zhu, PSERC ASU
%   and Ray Zimmerman, PSERC Cornell
%   Based on mpraw2mp.m, written by: Yujia Zhu, Jan 2014, yzhu54@asu.edu.
%   Copyright (c) 2014 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

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
data = psse_parse(records, sections, verbose, rev);

%% convert to MATPOWER case file
mpc = psse_convert(data);

%% (optionally) save MATPOWER case file
if ~isempty(mpc_name)
    if ~rev
        if data.id.REV
            rev = data.id.REV;
        else
            rev = 29;
        end
    end
    comments = { upper(mpc_name) };
    for k = 0:2
        str = data.id.(sprintf('comment%d', k));
        if ~isempty(str)
            comments{end+1} = sprintf('   %s', str);
        end
    end
    comments{end+1} = '';
    comments{end+1} = sprintf('   Converted by MATPOWER %s using PSSE2MPC on %s', mpver, date);
    comments{end+1} = sprintf('   from ''%s'' using PSS/E rev %d format.', rawfile_name, rev);

    if verbose
        spacers = repmat('.', 1, 45-length(mpc_name));
        fprintf('Saving to MATPOWER case ''%s'' %s', mpc_name, spacers);
    end
    savecase(mpc_name, comments, mpc);
    if verbose
        fprintf(' done.\n');
    end
end
