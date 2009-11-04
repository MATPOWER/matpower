function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
                runopf(casename, mpopt, fname, solvedcase)
%RUNOPF  Runs an optimal power flow.
%
%   Output arguments options:
%
%   results = runopf(...)
%   [results, success] = runopf(...)
%   [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(...)
%
%   Input arguments options:
%
%   runopf(casename)
%   runopf(casename, mpopt)
%   runopf(casename, mpopt, fname)
%   runopf(casename, mpopt, fname, solvedcase)
%
%   Runs an optimal power flow and optionally returns the solved values in
%   the data matrices, the objective function value, a flag which is true if
%   the algorithm was successful in finding a solution, and the elapsed time
%   in seconds. Alternatively, the solution can be returned as fields in a
%   results struct and an optional success flag.
%
%   All input arguments are optional. If casename is provided it specifies
%   the name of the input data file or struct (see also 'help caseformat' and
%   'help loadcase') containing the opf data. The default value is 'case9'. If
%   the mpopt is provided it overrides the default MATPOWER options vector and
%   can be used to specify the solution algorithm and output options among
%   other things (see 'help mpoption' for details). If the 3rd argument is
%   given the pretty printed output will be appended to the file whose name is
%   given in fname. If solvedcase is specified the solved case will be written
%   to a case file in MATPOWER format with the specified name. If solvedcase
%   ends with '.mat' it saves the case as a MAT-file otherwise it saves it as
%   an M-file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%-----  initialize  -----
%% default arguments
if nargin < 4
    solvedcase = '';                %% don't save solved case
    if nargin < 3
        fname = '';                 %% don't print results to a file
        if nargin < 2
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casename = 'case9'; %% default data file is 'case9.m'
            end
        end
    end
end

%%-----  run the optimal power flow  -----
[r, success] = opf(casename, mpopt);

%%-----  output results  -----
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        printpf(r, fd, mpopt);
        fclose(fd);
    end
end
printpf(r, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, r);
end

if nargout == 1 || nargout == 2
    MVAbase = r;
    bus = success;
elseif nargout > 2
    [MVAbase, bus, gen, gencost, branch, f, et] = ...
        deal(r.baseMVA, r.bus, r.gen, r.gencost, r.branch, r.f, r.et);
% else  %% don't define MVAbase, so it doesn't print anything
end
