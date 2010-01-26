function [MVAbase, bus, gen, gencost, branch, f, success, et] = ...
                runopf(casedata, mpopt, fname, solvedcase)
%RUNOPF  Runs an optimal power flow.
%
%   Runs an optimal power flow (AC OPF by default) optionally returning
%   a results struct and success flag.
%
%   results = runopf
%   results = runopf(casedata)
%   results = runopf(casedata, mpopt)
%   results = runopf(casedata, mpopt, fname)
%   results = runopf(casedata, mpopt, fname, solvedcase)
%   [results, success] = runopf(casedata, mpopt, fname, solvedcase)
%
%   Inputs (all are optional):
%       casedata : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also 'help caseformat' and 'help loadcase')
%       mpopt : MATPOWER options vector to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more.
%           (see also 'help mpoption')
%       fname : name of a file to which the pretty-printed output will
%           be appended
%       solvedcase : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       results : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = success, 0 = failure
%           (additional OPF fields, see 'help opf' for details)
%       success : the success flag can additionally be returned as
%           a second output argument
%
%   Alternatively, for compatibility with previous versions of MATPOWER,
%   some of the results can be returned as individual output arguments:
%       [baseMVA, bus, gen, gencost, branch, f, success, et] = runopf(...)

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

%%-----  run the optimal power flow  -----
[r, success] = opf(casedata, mpopt);

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
