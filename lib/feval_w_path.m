function varargout = feval_w_path(fpath, fname, varargin)
%FEVAL_W_PATH  Calls a function located by the specified path.
%   VARARGOUT = FEVAL_W_PATH(FPATH, FNAME, VARARGIN)
%
%   Calls a function whose path is given, even if the function is outside
%   of the Matlab path. Assumes that the current working directory is always
%   first in the Matlab path.
%
%   Inputs:
%       FPATH - string containing the path to the function to be called,
%               can be absolute or relative to current working directory
%       FNAME - string containing the name of the function to be called
%       VARARGIN - variable number of arguments to be passed to the function
%
%   Output:
%       VARARGOUT - variable number of return arguments (depends on the caller)
%
%   Note that any sub-functions located in the directory specified by FPATH
%   will also be available to be called by the FNAME function.
%
%   Examples:
%       % Assume '/opt/testfunctions' is NOT in the Matlab path.
%       rv = feval_w_path('/opt/testfunctions', 'mytestfcn', arg1, arg2);

%   MATPOWER
%   Copyright (c) 2016 by Power System Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   $Id$
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% check input type
if ~ischar(fpath)
    error('feval_w_path: FPATH must be a string');
end
if ~ischar(fname)
    error('feval_w_path: FNAME must be a string');
end

%% see if path exists
if exist(fpath, 'dir') ~= 7
    error('feval_w_path: Sorry, ''%s'' is not a valid directory path.', fpath);
end

cwd = pwd;      %% save the current working dir
cd(fpath);      %% switch to the dir with the mfile
try
    [varargout{1:nargout}] = feval(fname, varargin{:});
    cd(cwd);    %% switch back to saved dir
catch
    cd(cwd);    %% switch back to saved dir
    rethrow(lasterror);
end
