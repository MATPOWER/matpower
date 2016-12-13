function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
                runuopf(casedata, mpopt, fname, solvedcase)
%RUNUOPF  Runs an optimal power flow with unit-decommitment heuristic.
%   [RESULTS, SUCCESS] = RUNUOPF(CASEDATA, MPOPT, FNAME, SOLVEDCASE)
%
%   Runs an optimal power flow (AC OPF by default) with a heuristic which
%   allows it to shut down "expensive" generators, optionally returning
%   a RESULTS struct and SUCCESS flag.
%
%   Inputs (all are optional):
%       CASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
%       MPOPT : MATPOWER options struct to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more (see also MPOPTION).
%       FNAME : name of a file to which the pretty-printed output will
%           be appended
%       SOLVEDCASE : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       RESULTS : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = succeeded, 0 = failed
%           (additional OPF fields, see OPF for details)
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       results = runuopf;
%       results = runuopf(casedata);
%       results = runuopf(casedata, mpopt);
%       results = runuopf(casedata, mpopt, fname);
%       results = runuopf(casedata, mpopt, fname, solvedcase);
%       [results, success] = runuopf(...);
%
%       Alternatively, for compatibility with previous versions of MATPOWER,
%       some of the results can be returned as individual output arguments:
%
%       [baseMVA, bus, gen, gencost, branch, f, success, et] = runuopf(...);
%
%   Example:
%       results = runuopf('case30');
%
%   See also RUNOPF, RUNDUOPF.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
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
                casedata = 'case9'; %% default data file is 'case9.m'
            end
        end
    end
end

%%-----  run the unit de-commitment / optimal power flow  -----
[r, success] = uopf(casedata, mpopt);

%%-----  output results  -----
if fname
    [fd, msg] = fopen(fname, 'at');
    if fd == -1
        error(msg);
    else
        if mpopt.out.all == 0
            printpf(r, fd, mpoption(mpopt, 'out.all', -1));
        else
            printpf(r, fd, mpopt);
        end
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
