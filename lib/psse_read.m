function [records, sections] = psse_read(rawfile_name, verbose)
%PSSE_READ  Reads the data from a PSS/E RAW data file.
%   [RECORDS, SECTIONS] = PSSE_READ(RAWFILE_NAME)
%   [RECORDS, SECTIONS] = PSSE_READ(RAWFILE_NAME, VERBOSE)
%
%   Reads the data from a PSS/E RAW data file into a cell array of
%   strings, corresponding to the lines/records in the file. It
%   detects the beginning and ending indexes of each section as well
%   as any Q record used to indicate the end of the data.
%
%   Input:
%       RAWFILE_NAME :  name of the PSS/E RAW file to be read
%                       (opened directly with FILEREAD)
%       VERBOSE      :  1 to display progress info, 0 (default) otherwise
%
%   Output:
%       RECORDS :   a cell array of strings, one for each line in
%                   the file (new line characters not included)
%       SECTIONS  : a struct array with the following fields
%           first : index into RECORDS of first line of the section
%           last  : index into RECORDS of last line of the section
%           name  : name of the section (if available) extracted
%                   from the 'END OF <NAME> DATA, BEGIN ... DATA'
%                   comment typically found in the terminator line
%
%   See also PSSE2MPC.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
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

%% default args
if nargin < 2
    verbose = 0;
end

%% define some functions needed for use with CELLFUN
end_of_section  = @(s)~isempty(regexp(s, '^(Q|\s*0)[\s]*(/.*)?$', 'once'));
q_record        = @(s)~isempty(regexp(s, '^Q', 'once'));
section_name    = @(s)regexp(s, '/\s*END OF (.*)\s', 'tokens');

%% check for regexp split support
if ~have_fcn('regexp_split')
    error('psse_read: Sorry, but PSSE2MPC requires support for the ''split'' argument to regexp(), so it does not work on versions of Matlab prior to 7.3 or Octave prior to 3.8.');
end

%% read in the file, split on newlines
if verbose
    spacers = repmat('.', 1, 56-length(rawfile_name));
    fprintf('Reading file ''%s'' %s', rawfile_name, spacers);
end
str = fileread(rawfile_name);
if verbose
    fprintf(' done.\nSplitting into individual lines ...');
end
records = regexp(str, '\n|\r\n|\r', 'split');
if verbose
    str = sprintf('%d lines read', length(records));
    spacers = repmat('.', 1, 32-length(str));
    fprintf('%s %s ... done.\nAnalyzing sections ...', spacers, str);
end

%% find end of section and/or Q record
eos = find(cellfun(end_of_section, records));

%% terminate at Q record
qi = find(cellfun(q_record, records(eos)));
if ~isempty(qi)         %% remove everything beyond Q record ...
    if qi(1) > 1 && eos(qi(1)) - eos(qi(1)-1) > 1
        %% ... leaving Q record itself, if it is a section terminator
        records(eos(qi(1))+1:end) = [];
        eos(qi(1)+1:end) = [];
    else
        %% ... removing Q record too
        records(eos(qi(1)):end) = [];
        eos(qi(1):end) = [];
    end
end

%% find section extents and names
ns = length(eos) + 1;
i1 = [1 4 eos(1:end-1)+1];
iN = [3 eos-1];
names = cell(1, ns);
names{1} = 'ID';
for k = 2:ns
    tmp = regexp(records{eos(k-1)}, 'DATA', 'split');
    if isempty(tmp)     %% workaround a bug in Matlab 7.3 (R2006b) (and possibly earlier)
        names{k} = '';
    else
        tmp2 = section_name(tmp{1});
        if isempty(tmp2)
            names{k} = '';
        else
            names{k} = tmp2{1}{1};
        end
    end
end
if verbose
    str = sprintf('%d sections', ns);
    spacers = repmat('.', 1, 45-length(str));
    fprintf('%s %s ... done.\n', spacers, str);
end

%% create the sections struct
sections = struct(  'first', num2cell(i1), ...
                    'last', num2cell(iN), ...
                    'name', names   );
