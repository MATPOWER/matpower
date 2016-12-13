function [mpc, warnings] = cdf2mpc(cdf_file_name, mpc_name, verbose)
%CDF2MPC  Converts an IEEE CDF data file into a MATPOWER case struct.
%   MPC = CDF2MPC(CDF_FILE_NAME)
%   MPC = CDF2MPC(CDF_FILE_NAME, VERBOSE)
%   MPC = CDF2MPC(CDF_FILE_NAME, MPC_NAME)
%   MPC = CDF2MPC(CDF_FILE_NAME, MPC_NAME, VERBOSE)
%   [MPC, WARNINGS] = CDF2MPC(CDF_FILE_NAME, ...)
%
%   Converts an IEEE Common Data Format (CDF) data file into a MATPOWER case
%   struct.
%
%   Input:
%       CDF_FILE_NAME :  name of the IEEE CDF file to be converted
%       MPC_NAME      :  (optional) file name to use to save the resulting
%                         MATPOWER case
%       VERBOSE       :  1 (default) to display progress info, 0 otherwise
%
%   Output(s):
%       MPC      : resulting MATPOWER case struct
%       WARNINGS : (optional) cell array of strings containing warning
%                  messages (included by default in comments of MPC_NAME).
%
%   The IEEE CDF does not include some data need to run an optimal power
%   flow. This script creates default values for some of this data as
%   follows:
%
%       Bus data:
%           Vmin = 0.94 p.u.
%           Vmax = 1.06 p.u.
%       Gen data:
%           Pmin = 0 MW
%           Pmax = Pg + baseMVA
%       Gen cost data:
%           Quadratic costs with:
%               c2 = 10 / Pg, c1 = 20, c0 = 0, if Pg is non-zero, and
%               c2 = 0.01,    c1 = 40, c0 = 0, if Pg is zero
%           This should yield an OPF solution "close" to the
%           existing solution (assuming it is a solved case)
%           with lambdas near $40/MWh. See 'help caseformat'
%           for details on the cost curve format.
%
%   CDF2MPC may modify some of the data which are "infeasible" for
%   running optimal power flow. If so, warning information will be
%   printed out on screen.
%
%   Note: Since our code can not handle transformers with variable tap,
%   you may not expect to get exactly the same power flow solution
%   using converted data. This is the case when we converted ieee300.cdf.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Deqiang (David) Gan, PSERC Cornell & Zhejiang University
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% handle input args
if nargin < 2
    verbose = 1;
    mpc_name = '';
elseif ischar(mpc_name)     %% save the file
    if nargin < 3
        verbose = 1;
    end
else                        %% don't save the file
    verbose = mpc_name;
    mpc_name = '';
end

%%-----  read data from CDF file into mpc fields
%% verify valid input filename
[cdf_path cdf_name cdf_ext] = fileparts(cdf_file_name);
if isempty(cdf_ext)
    cdf_ext = '.cdf';
    cdf_file_name = strcat(cdf_name , cdf_ext);
end

%% open input file
[fid, msg] = fopen(cdf_file_name, 'r');
if fid < 0
    disp(msg);
    error('cdf2mpc: Can not read the input file:  %s', cdf_file_name )
end
if verbose
    fprintf('Converting file ''%s''\n', cdf_file_name);
    fprintf('  WARNINGS:\n');
end

%% initialize list of warnings
warnings = {};

% get baseMVA
title_cdf = fgetl(fid);
if isnumeric(str2num(title_cdf(32:37))) && not(isempty(str2num(title_cdf(32:37))))
    baseMVA = str2num(title_cdf(32:37));
    if length(findstr(title_cdf(2:9), '/')) == 2    %% date in the file 
        warnings{end+1} = sprintf('check the title format in the first line of the cdf file.');
    end
else
    error('cdf2mpc: Error getting the Base MVA, check the title format in the first line of the file.')
end
    
% find string 'BUS DATA FOLLOWS'
while 1
    line = fgetl(fid);
    if line(1:16) == 'BUS DATA FOLLOWS', break, end
end

% ----- get bus data, feed them into matrix bus, gen, gencost
ibus = 0;
igen = 0;
iarea = 0;

while 1
    line = fgetl(fid);
    if line(1:4) == '-999', break, end

    % feed bus data
    ibus = ibus + 1;
    bus(ibus, BUS_I) = str2num(line(1:4));          % bus number
    bus_name{ibus} = strtrim(line(6:17));           % bus names
    bus(ibus, BUS_TYPE) = str2num(line(25:26));
    if bus(ibus, BUS_TYPE) == 0                     % bus type
        bus(ibus, BUS_TYPE) = 1;
    end
    if (bus(ibus, BUS_TYPE) < 2)                    % Pd
        bus(ibus, PD) = str2num(line(41:49)) - str2num(line(60:67));
    elseif (bus(ibus, BUS_TYPE) >= 2)
        bus(ibus, PD) = str2num(line(41:49));
    end
    bus(ibus, QD) = str2num(line(50:59));           % Qd
    bus(ibus, GS) = baseMVA*str2num(line(107:114)); % Gs
    bus(ibus, BS) = baseMVA*str2num(line(115:122)); % Bs
    bus(ibus, BUS_AREA) = str2num(line(19:20));     % area
    bus(ibus, VM) = str2num(line(28:33));           % Vm
    bus(ibus, VA) = str2num(line(34:40));           % Va
    bus(ibus, BASE_KV) = str2num(line(77:83));      % baseKV
    bus(ibus, ZONE) = str2num(line(21:23));         % zone
    bus(ibus, VMAX) = 1.06;                         % default voltage upper limit
    bus(ibus, VMIN) = 0.94;                         % default voltage lower limit

    % feed gen and gencost
    Pg = str2num(line(60:67));
    Qg = str2num(line(68:75));
    Qmax = str2num(line(91:98));
    Qmin = str2num(line(99:106));
    if bus(ibus, BUS_TYPE) >= 2
        igen = igen + 1;
        if bus(ibus, BUS_TYPE) == 3, refgen = igen; end
        gen(igen, GEN_BUS) = bus(ibus, BUS_I);      % bus number
        gen(igen, PG) = Pg;                         % Pg
        if gen(igen, PG) < 0                        % negative Pg is transformed as load
            bus(ibus, PD) = bus(ibus, PD) - gen(igen, PG);
            warnings{end+1} = sprintf('negative Pg at bus %g treated as Pd', bus(ibus, BUS_I));
            if verbose
                fprintf('    %s\n', warnings{end});
            end
            gen(igen, PG) = 0;
        end
        gen(igen, QG)   = Qg;                       % Qg
        gen(igen, QMAX) = Qmax;                     % Qmax
        gen(igen, QMIN) = Qmin;                     % Qmin
        if Qmax - Qmin < 0.01                       % Qmax is modified
            gen(igen, QMAX) = Qmin + 0.1 * baseMVA;
            warnings{end+1} = sprintf('Qmax = Qmin at generator at bus %4i (Qmax set to Qmin + %g)', bus(ibus, BUS_I), baseMVA/10);
            if verbose
                fprintf('    %s\n', warnings{end});
            end
        end
        gen(igen, VG)    = str2num(line(85:90));    % specified voltage
        gen(igen, MBASE) = baseMVA;                 % baseMVA
        gen(igen, GEN_STATUS) = 1;                  % default status is 'on'
        gen(igen, PMAX)  = gen(igen, 2) + baseMVA;  % Pmax
        gen(igen, PMIN)  = 0;                       % Pmin = 0 by default

        gencost(igen, MODEL)    = POLYNOMIAL;       % by default, sets the model as polynomial
        gencost(igen, STARTUP)  = 0;                % start up cost is zero by default
        gencost(igen, SHUTDOWN) = 0;                % shut down cost is zero by default
        gencost(igen, NCOST)    = 3;                % number of coefficients in polynomial cost
%       gencost(igen, COST)     = 0.01;             % default c2
%       gencost(igen, COST+1)   = 0.3;              % default c1
%       gencost(igen, COST+2)   = 0.2;              % default c0
    end
end

totload = sum(bus(:, PD));
totgen = sum(gen(:, PG));
if totgen < 1.04 * totload
    gen(refgen, PMAX) = gen(refgen, PG) + 1.1 * totload - totgen;   % Pg at slack bus is modified
    warnings{end+1} = sprintf('Insufficient generation, setting Pmax at slack bus (bus %d) to %g', gen(refgen, [GEN_BUS, PMAX]));
    if verbose
        fprintf('    %s\n', warnings{end});
    end
end

% ----- set up the cost coefficients of generators
ng = size(gen, 1);
% gencost(:, COST)    = zeros(ng, 1);
% gencost(:, COST+1)  = 100*ones(ng, 1) ./ (gen(:, PG) + 10*ones(ng, 1));
% gencost(:, COST+2)  = 100*ones(ng, 1) ./ (gen(:, PG) + 10*ones(ng, 1));
zg  = find(gen(:, PG) == 0);                %% for Pg = 0
gencost(zg, COST)  = 0.01 * ones(size(zg));
gencost(zg, COST+1) = 40 * ones(size(zg));
nzg = find(gen(:, PG) ~= 0);                %% Pg non-zero
gencost(nzg, COST) = 10 * ones(size(nzg)) ./ gen(nzg, PG);
gencost(nzg, COST+1) = 20 * ones(size(nzg));
gencost(:, COST+2) = zeros(ng, 1);

% find string 'BRANCH DATA FOLLOWS'
while 1
    line = fgetl(fid);
    if line(1:19) == 'BRANCH DATA FOLLOWS', break, end
end

% ----- get branch data, feed them into matrix branch
k = 0;
while 1
    line = fgetl(fid);
    if line(1:4) == '-999', break, end

    k = k + 1;
    branch(k, F_BUS)  = str2num(line(1:4));     % fbus (also the tap bus)
    branch(k, T_BUS)  = str2num(line(6:9));     % tbus
    branch(k, BR_R)   = str2num(line(20:29));   % R
    branch(k, BR_X)   = str2num(line(30:40));   % X
    branch(k, BR_B)   = str2num(line(41:50));   % B
    branch(k, RATE_A) = str2num(line(51:55));   % RATE A
    if branch(k, RATE_A) < 0.000001
        branch(k, RATE_A) = 0;                  % RATE A is modified
        warnings{end+1} = sprintf('MVA limit of branch %d - %d not given, set to %g', branch(k, [F_BUS, T_BUS, RATE_A]));
        if verbose
            fprintf('    %s\n', warnings{end});
        end
    end
    branch(k, RATE_B) = str2num(line(57:61));   % RATE B
    branch(k, RATE_C) = str2num(line(63:67));   % RATE C
    branch(k, TAP)    = str2num(line(77:82));   % transformer turns ratio
    branch(k, SHIFT)  = 0;                      % phase shifter can not be modelled
    branch(k, BR_STATUS) = 1;                   % by default, branch is on
end
if verbose
    fprintf('Done.\n');
end
fclose(fid);

% put in struct
mpc.baseMVA = baseMVA;
mpc.bus     = bus;
mpc.branch  = branch;
mpc.gen     = gen;
mpc.gencost = gencost;
mpc.bus_name = bus_name;
mpc = loadcase(mpc);    %% convert to internal (e.g. v. '2') case format

%% (optionally) save MATPOWER case file
if ~isempty(mpc_name)
    comments = {''};
    if ~isempty(title_cdf)
        comments{end+1} = sprintf('   %s', title_cdf);
    end
    comments{end+1} = '';
    comments{end+1} = sprintf('   Converted by MATPOWER %s using CDF2MPC on %s', mpver, date);
    comments{end+1} = sprintf('   from ''%s''.', cdf_file_name);

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
