function varargout = feval_w_path(fpath, fname, varargin)
%FEVAL_W_PATH  Calls a function located by the specified path.
%   FEVAL_W_PATH(FPATH, F, x1, ..., xn)
%   [y1, ..., yn] = FEVAL_W_PATH(FPATH, F, x1, ..., xn)
%
%   Identical to the built-in FEVAL, except that the function F need not be
%   in the MATLAB/Octave path if it is defined in a file in the path specified
%   by FPATH. Assumes that the current working directory is always first in
%   the MATLAB/Octave path.
%
%   Inputs:
%       FPATH - string containing the path to the function to be called,
%               can be absolute or relative to current working directory
%       F - string containing the name of the function to be called
%       x1, ..., xn - variable number of input arguments to be passed to F
%
%   Output:
%       y1, ..., yn - variable number arguments returned by F (depending on
%                     the caller)
%
%   Note that any sub-functions located in the directory specified by FPATH
%   will also be available to be called by the F function.
%
%   Examples:
%       % Assume '/opt/testfunctions' is NOT in the MATLAB/Octave path, but
%       % /opt/testfunctions/mytestfcn.m defines the function mytestfcn()
%       % which takes 2 input arguments and outputs 1 return argument.
%       y = feval_w_path('/opt/testfunctions', 'mytestfcn', x1, x2);

%   MATPOWER
%   Copyright (c) 2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
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

if isempty(fpath)   %% just call feval directly, no cd necessary
    [varargout{1:nargout}] = feval(fname, varargin{:});
else
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
end
