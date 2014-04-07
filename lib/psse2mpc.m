function mpc = psse2mpc(rawfile_name)
%PSSE2MPC  Converts a PSS/E RAW data file into a MATPOWER case struct.
%   MPC = PSSE2MPC(RAWFILE_NAME)
%
%   Converts a PSS/E RAW data file into a MATPOWER case struct.
%
%   Input:
%       RAWFILE_NAME : the name of the PSS/E RAW file to be converted
%           (opened directly with fopen)
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

%% read data from version 33 PSS/E RAW file
data = psse_read_33(rawfile_name);

%% convert to MATPOWER case file
mpc = psse_convert_33(data);
