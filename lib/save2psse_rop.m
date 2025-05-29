% This is a new function (not part of MATPOWER) developed to generate a 
% PSS/E ROP (Raw Operating Point) file from a given test case.
function fname_out = save2psse_rop(fname, mpc, rawver)
%SAVE2PSSE_ROP   Save MATPOWER case to a PSS/E ROP (Raw Operating Point) file.
%
%   FNAME_OUT = SAVE2PSSE_ROP(FNAME, MPC, ROPVER)
%
%   Inputs:
%       FNAME   - Filename (string) to save the ROP file to. If no extension is
%                 provided, '.rop' is appended automatically.
%       MPC     - MATPOWER case struct containing bus, gen, branch, gencost matrices.
%       ROPVER  - Version number for the PSS/E ROP format (currently not used).
%
%   Output:
%       FNAME_OUT - Full path of the output file (with extension).
%
%   This function exports the data from a MATPOWER case into a text file
%   formatted in the ROP style compatible with PSS/E. Currently, it includes:
%       - Generator dispatch data
%       - Active power dispatch tables
%       - Polynomial cost tables
%       - Section separators for other categories (left empty)
%
%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   and individual contributors (see AUTHORS file for details).
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, C0, C1, C2] = idx_quad_cost;
c = idx_dcline;

%% verify valid filename
[pathstr, fcn_name, extension] = fileparts(fname);
if isempty(extension)
    extension = '.raw';
end
fname = fullfile(pathstr, [fcn_name extension]);

%% open file
[fd, msg] = fopen(fname, 'wt');     %% print it to a text file
if fd == -1
    error(['save2psse: ', msg]);
end

%% create map of external bus numbers to bus indices
i2e = mpc.bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(mpc.bus, 1))';

%% Case Identification Data
fprintf(fd, ' 0          /  PSS/E-30.0, %s\n', datestr(clock));
%% Bus Voltage Constraint data
fprintf(fd, ' 0 / End of Bus Voltage Constraint data, begin Adjustable Bus Shunt data\n');
%% Adjustable Bus Shunt data
fprintf(fd, ' 0 / End of Adjustable Bus Shunt data, begin Bus Load data\n');
%% Bus Load data
fprintf(fd, ' 0 / End of Bus Load data, begin Adjustable Bus Load Tables\n');
%% Adjustable Bus Load Tables
fprintf(fd, ' 0 / End of Adjustable Bus Load Tables, begin Generator Dispatch data\n');

%% Gen Data
ng = size(mpc.gen, 1);  %% number of buses
%% Generator Dispatch data
%% Bus, GenID, Disp, DspTbl
GenIDs = zeros(length(mpc.bus),1);
for i = 1:ng
    GenIDs(mpc.gen(i, GEN_BUS)) = GenIDs(mpc.gen(i, GEN_BUS)) + 1;
    GenID = GenIDs(mpc.gen(i, GEN_BUS));
    fprintf(fd, '%d, %d, %.1f, %d\n', ...
    [ mpc.gen(i, GEN_BUS), ...
    GenID, ...
    1, ...
    i ...
    ]');
end

fprintf(fd, ' 0 / End of Generator Dispatch data, begin Active Power Dispatch Tables\n');

%% Active Power Dispatch Tables
%% Tbl, P_max, P_min, Fuel_cost, C_typ, Status, C_tbl
fprintf(fd, '%d, %.6f, %.6f, %.1f, %d, %d, %d\n', ...
    [ (1:ng)', ...
    mpc.gen(:, [PMAX]), ...
    mpc.gen(:, [PMIN]), ...
    ones(ng, 1), ...
    ones(ng, 1), ...
    ones(ng, 1), ...
    (1:ng)' ...
    ]');
fprintf(fd, ' 0 / End of Active Power Dispatch Tables, begin Generation Reserve data\n');
%% Generation Reserve data
fprintf(fd, ' 0 / End of Generation Reserve data, begin Generation Reactive Capability data\n');
%% Generation Reactive Capability data
fprintf(fd, ' 0 / End of Generation Reactive Capability data, begin Adjustable Branch Reactance data\n');
%% Adjustable Branch Reactance data
fprintf(fd, ' 0 / End of Adjustable Branch Reactance data, begin Piece - wise Linear Cost Tables\n');
%% Piece-wise Linear Cost Tables
fprintf(fd, ' 0 / End of Piece-wise Linear Cost Tables, begin Piece-wise Quadratic Cost Tables\n');
%% Piece-wise Quadratic Cost Tables
fprintf(fd, ' 0 / End of Piece-wise Quadratic Cost Tables, begin Polynomial Cost Tables\n');
%% Polynomial Cost Tables
%%  Ptbl, Label, Cost, Costlin, Costquad, Costexp, Expn
gl = cell2mat(cellfun(@(s) sprintf('QUADRATIC %-4d', s), ...
    num2cell((1:ng)'), 'UniformOutput', 0));
if mpc.gencost(1,NCOST) == 3
    for i = 1:ng
        % fprintf(fd, '%d, %s, %.6f, %.6f, %.6f, %.6f, %.6f\n', ...
        fprintf(fd, '%d, ''%c%c%c%c%c%c%c%c%c%c%c%c%c%c'', %.6f, %.6f, %.6f, %.1f, %.1f\n', ...
            i, ...
            gl(i, :), ...
            mpc.gencost(i, C0), ...
            mpc.gencost(i, C1), ...
            mpc.gencost(i, C2), ...
            0.0, ...
            0.0);
    end
else
    for i = 1:ng
        % fprintf(fd, '%d, %s, %.6f, %.6f, %.6f, %.6f, %.6f\n', ...
        fprintf(fd, '%d, ''%c%c%c%c%c%c%c%c%c%c%c%c%c%c'', %.1f, %.6f, %.6f, %.1f, %.1f\n', ...
            i, ...
            gl(i, :), ...
            0.0, ...
            mpc.gencost(i, C1), ...
            mpc.gencost(i, C2), ...
            0.0, ...
            0.0);
    end
end
fprintf(fd, ' 0 / End of Polynomial Cost Tables, begin Period Reserve data\n');
%% Period Reserve data
fprintf(fd, ' 0 / End of Period Reserve data, begin Branch Flow Constraint data\n');
%% Branch Flow Constraint data
fprintf(fd, ' 0 / End of Branch Flow Constraint data, begin Interface Flow data\n');
%% Interface Flow data
fprintf(fd, ' 0 / End of Interface Flow data, begin Linear Constraint Equation Dependency data\n');
%% Linear Constraint Equation Dependency data
fprintf(fd, ' 0 / End of Linear Constraint Equation Dependency data, begin 2-terminal dc Line Constraint data\n');
%% 2-terminal dc Line Constraint data
fprintf(fd, ' 0 / End of 2-terminal dc Line Constraint data\n');

%% close file
if fd ~= 1
    fclose(fd);
end

if nargout > 0
    fname_out = fname;
end

function ckt = generate_ckt_num(ft, ckt)
% returns a vector of circuit numbers given a 2-col matrix of from/to bus nums
% where duplicate pairs get incremented ckt numbers
n = size(ft, 1);
if nargin < 2
    ckt = zeros(n, 1);
end
ckt = ckt + 1;
[Au, iA, iAu] = unique(ft, 'rows', 'first');
k = find(~ismember([1:n]', iA));
if k
    ckt(k) = generate_ckt_num(ft(k, :), ckt(k));
end

function [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, C0, C1, C2] = idx_quad_cost
%IDX_COST   Defines constants for named column indices to gencost matrix.


%% define cost models
PW_LINEAR   = 1;
POLYNOMIAL  = 2;

%% define the indices
MODEL       = 1;    %% cost model, 1 = piecewise linear, 2 = polynomial
STARTUP     = 2;    %% startup cost in US dollars
SHUTDOWN    = 3;    %% shutdown cost in US dollars
NCOST       = 4;    %% number N = n+1 of end/breakpoints in piecewise linear
%% cost function, or of coefficients in polynomial cost fcn
COST        = 5;    %% parameters defining total cost function begin in this col
C2          = 5;                    %% (MODEL = 2) : cn, ..., c1, c0
C1          = 6;          %%      N coefficients of an n-th order polynomial cost fcn,
C0          = 7;         %%      starting with highest order, where cost is
%%      f(p) = cn*p^n + ... + c1*p + c0

