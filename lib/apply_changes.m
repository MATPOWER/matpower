function mpc = apply_changes(label, mpc, chgtab)
%APPLY_CHANGES  Applies a set of changes to a MATPOWER case
%   mpc_modified = apply_changes(label, mpc_original, chgtab)
%
%   Applies the set of changes identified by LABEL to the case in MPC, where
%   the change sets are specified in CHGTAB.
%
%   LABEL is an integer which identifies the set of changes of interest
%
%   MPC is a MATPOWER case struct with at least fields bus, gen and branch
%
%   CHGTAB is the table of changes that lists the individual changes required
%   by each change set, one "change" per row (multiple rows or changes
%   allowed for each change set). Type "help idx_ct" for more complete
%   information about the format.
%   1st column: change set label, integer > 0
%   2nd column: change set probability
%   3rd column: table to be modified (1:bus, 2:gen, 3:branch) or
%               4: bus area changes, apply to all gens/buses/branches in
%               a given area; 5: gen table area changes apply to all
%               generators in a given area; or 6: branch area changes apply
%               to all branches connected to buses in a given area.
%   4th column: row of table to be modified (if 3rd col is 1-3),
%               or area number for overall changes (if 3rd column is 4-6).
%   5th column: column of table to be modified. It is best to use the
%               named column index defined in the corresponding idx_bus,
%               idx_gen, idx_branch or idx_cost M-files in MATPOWER.
%   6th column: type of change: 1: absolute (replace)
%                               2: relative (multiply by factor)
%                               3: additive (add to value)
%   7th column: new value or multiplicative or additive factor
%
%   Examples:
%
%   chgtab = [ ...
%       1   0.1   CT_TGEN       2  GEN_STATUS     CT_REP  0;
%       2   0.05  CT_TGEN       3  PMAX           CT_REP  100;
%       3   0.2   CT_TBRCH      2  BR_STATUS      CT_REP  0;
%       4   0.1   CT_TAREALOAD  2  CT_LOAD_ALL_P  CT_REL  1.1;
%    ];
%
%   Description of each change set:
%   1. Turn off generator 2, 10% probability.
%   2. Set generator 3's max output to 100 MW, 5% probability.
%   3. Take branch 2 out of service, 20% probability.
%   4. Scale all loads in area 2 (real & reactive, fixed and dispatchable)
%      by a factor of 1.1, 10% probability.

%   To do:
%       - check for valid row number

%   MATPOWER
%   Copyright (c) 2000-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% find the change set to be applied
kk = find(label == chgtab(:, CT_LABEL));
if isempty(kk)
    error('apply_changes: LABEL %d not found in CHGTAB', label);
end

%% create map of external bus numbers to bus indices
i2e = mpc.bus(:, BUS_I);
e2i = sparse(max(i2e), 1);
e2i(i2e) = (1:size(mpc.bus, 1))';

for c = 1:length(kk)
    %% extract data for individual change to be applied
    %% (c-th change in specified change set)
    change = num2cell(chgtab(kk(c), :));
    [tbl, row, col, typ, val] = ...
        deal(change{[CT_TABLE, CT_ROW, CT_COL, CT_CHGTYPE, CT_NEWVAL]});

    %% apply the change
    if tbl == CT_TBUS                               %% modify bus table
      if ~any(col == [PD QD GS BS VMAX VMIN])
        error('apply_changes: modification to column %d of bus table not supported', col);
      end
      if row == 0                                       %% modify all rows
        if typ == CT_REP                                    %% replace
          mpc.bus(:, col) = val * ones(size(mpc.bus, 1), 1);
        elseif typ == CT_REL                                %% scale
          mpc.bus(:, col) = val * mpc.bus(:, col);
        elseif typ == CT_ADD                                %% shift
          mpc.bus(:, col) = val + mpc.bus(:, col);
        else
          error('apply_changes: unsupported modification type %d for bus table', typ);
        end
      else                                              %% modify single row
        if typ == CT_REP                                    %% replace
          mpc.bus(row, col) = val;
        elseif typ == CT_REL                                %% scale
          mpc.bus(row, col) = val * mpc.bus(row, col);
        elseif typ == CT_ADD                                %% shift
          mpc.bus(row, col) = val + mpc.bus(row, col);             
        else
          error('apply_changes: unsupported modification type %d for bus table', typ);
        end
      end
    elseif tbl == CT_TBRCH                          %% modify branch table
      if ~any(col == ...
            [BR_R BR_X BR_B RATE_A RATE_B RATE_C TAP SHIFT BR_STATUS ANGMIN ANGMAX])
        error('apply_changes: modification to column %d of branch table not supported', col);
      end
      if row == 0                                       %% modify all rows
        if typ == CT_REP                                    %% replace
          mpc.branch(:, col) = val * ones(size(mpc.branch, 1), 1);
        elseif typ == CT_REL                                %% scale
          mpc.branch(:, col) = val * mpc.branch(:, col);
        elseif typ == CT_ADD                                %% shift
          mpc.branch(:, col) = val + mpc.branch(:, col);
        else
          error('apply_changes: unsupported modification type %d for branch table', typ);
        end
      else                                              %% modify single row
        if typ == CT_REP                                    %% replace
          mpc.branch(row, col) = val;
        elseif typ == CT_REL                                %% scale
          mpc.branch(row, col) = val * mpc.branch(row, col);
        elseif typ == CT_ADD                                %% shift
          mpc.branch(row, col) = val + mpc.branch(row, col);
        else
          error('apply_changes: unsupported modification type %d for branch table', typ);
        end
      end
    elseif tbl == CT_TGEN                           %% modify gen table
      if ~any(col == ...
            [QMAX QMIN GEN_STATUS PMAX PMIN PC1 PC2 QC1MIN QC1MAX QC2MIN  QC2MAX ...
            RAMP_AGC RAMP_10 RAMP_30 RAMP_Q APF])
        error('apply_changes: modification to column %d of gen table not supported', col);
      end
      if row == 0                                       %% modify all rows
        if typ == CT_REP                                    %% replace
          mpc.gen(:, col) = val * ones(size(mpc.gen, 1), 1);
        elseif typ == CT_REL                                %% scale
          mpc.gen(:, col) = val * mpc.gen(:, col);
        elseif typ == CT_ADD                                %% shift
          mpc.gen(:, col) = val + mpc.gen(:, col);
        else
          error('apply_changes: unsupported modification type %d for gen table', typ);
        end
      else                                              %% modify single row
        if typ == CT_REP                                    %% replace
          mpc.gen(row, col) = val;
        elseif typ == CT_REL                                %% scale
          mpc.gen(row, col) = val * mpc.gen(row, col);
        elseif typ == CT_ADD                                %% shift
          mpc.gen(row, col) = val + mpc.gen(row, col);
        else
          error('apply_changes: unsupported modification type %d for gen table', typ);
        end
      end
    elseif tbl == CT_TGENCOST                       %% modify gencost table
      if col == CT_MODCOST_F || col == CT_MODCOST_X     %% use modcost to scale/shift cost
        if typ == CT_REL
          if col == CT_MODCOST_F
            modcost_type = 'SCALE_F';
          else  %% col == CT_MODCOST_X
            modcost_type = 'SCALE_X';
          end
        elseif typ == CT_ADD
          if col == CT_MODCOST_F
            modcost_type = 'SHIFT_F';
          else  %% col == CT_MODCOST_X
            modcost_type = 'SHIFT_X';
          end
        else
            error('apply_changes: unsupported modification type %d for gencost table CT_MODCOST_F/X modification', typ);
        end
        if row == 0
          mpc.gencost = modcost(mpc.gencost, val, modcost_type);
        else
          mpc.gencost(row, :) = modcost(mpc.gencost(row, :), val, modcost_type);
        end
      else                                              %% normal individual column mod
        if col < 1 || fix(col) ~= col   %% needs to be positive integer
          error('apply_changes: modification to column %d of gencost table not supported', col);
        end
        if row == 0                                       %% modify all rows
          if typ == CT_REP                                    %% replace
            mpc.gencost(:, col) = val * ones(size(mpc.gencost, 1), 1);
          elseif typ == CT_REL                                %% scale
            mpc.gencost(:, col) = val * mpc.gencost(:, col);
          elseif typ == CT_ADD                                %% shift
            mpc.gencost(:, col) = val + mpc.gencost(:, col);
          else
            error('apply_changes: unsupported modification type %d for gencost table', typ);
          end
        else                                              %% modify single row
          if typ == CT_REP                                    %% replace
            mpc.gencost(row, col) = val;
          elseif typ == CT_REL                                %% scale
            mpc.gencost(row, col) = val * mpc.gencost(row, col);
          elseif typ == CT_ADD                                %% shift
            mpc.gencost(row, col) = val + mpc.gencost(row, col);
          else
            error('apply_changes: unsupported modification type %d for gencost table', typ);
          end
        end
      end
    elseif tbl == CT_TAREABUS                   %% area-wide mod to bus table
      if ~any(col == [PD QD GS BS VMAX VMIN])
        error('apply_changes: area-wide modification to column %d of bus table not supported', col);
      end
      jj = find(mpc.bus(:, BUS_AREA) == row);
      if typ == CT_REP                                      %% replace
        mpc.bus(jj, col) = val * ones(size(jj));
      elseif typ == CT_REL                                  %% scale
        mpc.bus(jj, col) = val * mpc.bus(jj, col);
      elseif typ == CT_ADD                                  %% shift
        mpc.bus(jj, col) = val + mpc.bus(jj, col);
      else
        error('apply_changes: unsupported area-wide modification type %d for bus table', typ);
      end
    elseif tbl == CT_TAREABRCH                  %% area-wide mod to branch table
      if ~any(col == ...
            [BR_R BR_X BR_B RATE_A RATE_B RATE_C TAP SHIFT BR_STATUS ANGMIN ANGMAX])
        error('apply_changes: area-wide modification to column %d of branch table not supported', col);
      end
      jj = find((row == mpc.bus(e2i(mpc.branch(:, F_BUS)), BUS_AREA)) | ...
                (row == mpc.bus(e2i(mpc.branch(:, T_BUS)), BUS_AREA)) );
      if typ == CT_REP                                      %% replace
        mpc.branch(jj, col) = val * ones(size(jj));
      elseif typ == CT_REL                                  %% scale
        mpc.branch(jj, col) = val * mpc.branch(jj, col);
      elseif typ == CT_ADD                                  %% shift
        mpc.branch(jj, col) = val + mpc.branch(jj, col);
      else
        error('apply_changes: unsupported area-wide modification type %d for branch table', typ);
      end
    elseif tbl == CT_TAREAGEN                   %% area-wide mod to gen table
      if ~any(col == ...
            [QMAX QMIN GEN_STATUS PMAX PMIN PC1 PC2 QC1MIN QC1MAX QC2MIN  QC2MAX ...
            RAMP_AGC RAMP_10 RAMP_30 RAMP_Q APF])
        error('apply_changes: area-wide modification to column %d of gen table not supported', col);
      end
      jj = find(row == mpc.bus(e2i(mpc.gen(:, GEN_BUS)), BUS_AREA));
      if typ == CT_REP                                      %% replace
        mpc.gen(jj, col) = val * ones(size(jj));
      elseif typ == CT_REL                                  %% scale
        mpc.gen(jj, col) = val * mpc.gen(jj, col);
      elseif typ == CT_ADD                                  %% shift
        mpc.gen(jj, col) = val + mpc.gen(jj, col);
      else
        error('apply_changes: unsupported area-wide modification type %d for gen table', typ);
      end
    elseif tbl == CT_TAREAGENCOST               %% area-wide mod to gencost table
      jj = find(row == mpc.bus(e2i(mpc.gen(:, GEN_BUS)), BUS_AREA));
      if col == CT_MODCOST_F || col == CT_MODCOST_X
        if typ == CT_REL
          if col == CT_MODCOST_F
            modcost_type = 'SCALE_F';
          else  %% col == CT_MODCOST_X
            modcost_type = 'SCALE_X';
          end
        elseif typ == CT_ADD
          if col == CT_MODCOST_F
            modcost_type = 'SHIFT_F';
          else  %% col == CT_MODCOST_X
            modcost_type = 'SHIFT_X';
          end
        else
            error('apply_changes: unsupported area-wide modification type %d for gencost table CT_MODCOST_F/X modification', typ);
        end
        mpc.gencost(jj, :) = modcost(mpc.gencost(jj, :), val, modcost_type);
      else
        if col < 1 || fix(col) ~= col   %% needs to be positive integer
          error('apply_changes: area-wide modification to column %d of gencost table not supported', col);
        end
        if typ == CT_REP                                      %% replace
          mpc.gencost(jj, col) = val * ones(size(jj));
        elseif typ == CT_REL                                  %% scale
          mpc.gencost(jj, col) = val * mpc.gencost(jj, col);
        elseif typ == CT_ADD                                  %% shift
          mpc.gencost(jj, col) = val + mpc.gencost(jj, col);
        else
          error('apply_changes: unsupported area-wide modification type %d for gencost table', typ);
        end
      end
    elseif tbl == CT_TLOAD || tbl == CT_TAREALOAD           %% modify loads
      if ~any(abs(col) == (CT_LOAD_ALL_PQ:CT_LOAD_DIS_P))
        error('apply_changes: column=%d for load modifications is not supported', col);
      end
      switch abs(col)
        case {CT_LOAD_ALL_PQ, CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ}
            opt.pq = 'PQ';
        case {CT_LOAD_ALL_P, CT_LOAD_FIX_P, CT_LOAD_DIS_P}
            opt.pq = 'P';
      end
      switch abs(col)
        case {CT_LOAD_ALL_PQ, CT_LOAD_ALL_P}
            opt.which = 'BOTH';
        case {CT_LOAD_FIX_PQ, CT_LOAD_FIX_P}
            opt.which = 'FIXED';
        case {CT_LOAD_DIS_PQ, CT_LOAD_DIS_P}
            opt.which = 'DISPATCHABLE';
      end

      %% define load_zone
      if tbl == CT_TLOAD
        nb = size(mpc.bus, 1);
        if row == 0                                     %% modify load at all buses
          load_zone = ones(nb, 1);
        else                                            %% modify load at single bus
          load_zone = zeros(nb, 1);
          load_zone(row) = 1;
        end
      elseif tbl == CT_TAREALOAD
        load_zone = double(mpc.bus(:, BUS_AREA) == row);  %% modify load in single area
      end

      switch typ
        case CT_REP                                     %% replace
            opt.scale = 'QUANTITY';
            dmd = val;
        case CT_REL                                     %% scale
            opt.scale = 'FACTOR';
            dmd = val;
        case CT_ADD                                     %% shift
            opt.scale = 'QUANTITY';
            old_val = total_load(mpc, load_zone);
            dmd = old_val + val;
        otherwise
          error('apply_changes: unsupported modification type %d for loads', typ);
      end
      %% scale the loads ...
      if col < 0    %% ... including dispatchable load costs
          [mpc.bus, mpc.gen, mpc.gencost] = scale_load(dmd, mpc.bus, mpc.gen, load_zone, opt, mpc.gencost);
      else          %% ... not including dispatchable load costs
          [mpc.bus, mpc.gen] = scale_load(dmd, mpc.bus, mpc.gen, load_zone, opt);
      end
    else
      error('apply_changes: CHGTAB attempts to modify unsupported table type')
    end
end
