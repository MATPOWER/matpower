function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
                runuopf(casename, mpopt, fname, solvedcase)
%RUNUOPF  Runs an optimal power flow with unit-decommitment heuristic.
%
%   [baseMVA, bus, gen, gencost, branch, f, success, et] = ...
%           runuopf(casename, mpopt, fname, solvedcase)
%
%   Runs an optimal power flow with a heuristic which allows it to shut down
%   "expensive" generators and optionally returns the solved values in the
%   data matrices, the objective function value, a flag which is true if the
%   algorithm was successful in finding a solution, and the elapsed time in
%   seconds. All input arguments are optional. If casename is provided it
%   specifies the name of the input data file or struct (see also 'help
%   caseformat' and 'help loadcase') containing the opf data. The default
%   value is 'case9'. If the mpopt is provided it overrides the default
%   MATPOWER options vector and can be used to specify the solution
%   algorithm and output options among other things (see 'help mpoption' for
%   details). If the 3rd argument is given the pretty printed output will be
%   appended to the file whose name is given in fname. If solvedcase is
%   specified the solved case will be written to a case file in MATPOWER
%   format with the specified name. If solvedcase ends with '.mat' it saves
%   the case as a MAT-file otherwise it saves it as an M-file.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
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

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch, areas, gencost] = loadcase(casename);
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

%% run unit commitment / optimal power flow
[bus, gen, branch, f, success, et] = uopf(baseMVA, bus, gen, gencost, branch, areas, mpopt);

%% convert back to original bus numbering & print results
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        printpf(baseMVA, bus, gen, branch, f, success, et, fd, mpopt);
        fclose(fd);
    end
end
printpf(baseMVA, bus, gen, branch, f, success, et, 1, mpopt);

%% save solved case
if solvedcase
    savecase(solvedcase, baseMVA, bus, gen, branch, areas, gencost);
end

%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
if nargout, MVAbase = baseMVA; end

return;
