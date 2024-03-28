function t_most_tlmp(quiet)
% t_most_tlmp - Tests of TLMP for ramping and storage

%   MOST
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    quiet = 0;      %% verbose by default
end

t_begin(40, quiet);

if quiet
    verbose = 0;
else
    verbose = 0;
end
verbose = 0;

%% set up options
define_constants;
mpopt = mpoption('verbose', verbose);
mpopt = mpoption(mpopt, 'out.gen', 1);
mpopt = mpoption(mpopt, 'model', 'DC');
% mpopt = mpoption(mpopt, 'most.solver', 'GLPK');
% mpopt = mpoption(mpopt, 'most.solver', 'MOSEK');
mpopt = mpoption(mpopt, 'most.dc_model', 1);
if have_feature('mosek')
    sc = mosek_symbcon;
    mpopt = mpoption(mpopt, 'mosek.lp_alg', sc.MSK_OPTIMIZER_DUAL_SIMPLEX);     %% dual simplex
end
if ~verbose
    mpopt = mpoption(mpopt, 'out.all', 0);
end

%%-----  TLMP Example 1  -----
t = 'example 1 : ';
if verbose > 1
    fprintf('\n--------------------  TLMP Example 1  --------------------\n');
end
mpc = loadcase('t_case_tlmp');
xgd_table.colnames = {
    'PositiveLoadFollowReservePrice', ...
            'PositiveLoadFollowReserveQuantity', ...
                'NegativeLoadFollowReservePrice', ...
                        'NegativeLoadFollowReserveQuantity', ...
};
xgd_table.data = [
    1e-6    500 1e-6    500;
    1e-6    50  1e-6    50;
    1e-6    Inf 1e-6    Inf;
];
load_profile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REP, ...
    'values', [] );
load_profile.values(:, 1, 1) = [
    420;
    590;
    590;
]-1e-3;
ePg = [ 380 500 500; [40 90 90]-1e-3; -load_profile.values(:, 1, 1)' ];
eLMP = [25 35 30; 25 35 30; 25 35 30];
eTLMP = [25 35 30; 30 30 30; 25 35 30];

xgd = loadxgendata(xgd_table, mpc);
profiles = getprofiles(load_profile);
nt = size(profiles(1).values, 1);       % number of periods
mdi = loadmd(mpc, nt, xgd, [], [], profiles);
mdo = most(mdi, mpopt);
ms = most_summary(mdo);

t_is(mdo.results.success, 1, 12, [t 'success']);
t_is(ms.Pg, ePg, 12, [t 'Pg']);
t_is(mdo.results.GenPrices, eLMP, 5, [t 'LMP']);
t_is(mdo.results.GenTLMP, eTLMP, 5, [t 'TLMP']);

if verbose
    fprintf('    ');
    for t = 1:nt
        fprintf('%-20s', sprintf('   T=%d', t));
    end
    fprintf('\n');
    for g = 1:2
        fprintf('G%d  ', g);
        for t = 1:nt
            fprintf('%-20s', ...
                sprintf('(%.4g, %.4g, %.4g)', ms.Pg(g, t), ...
                                        mdo.results.GenPrices(g, t), ...
                                        mdo.results.GenTLMP(g, t) ));
        end
        fprintf('\n');
    end
end

%%-----  TLMP Example 2  -----
t = 'example 2';
if verbose > 1
    fprintf('\n--------------------  TLMP Example 2  --------------------\n');
end
load_profile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REP, ...
    'values', [] );
load_profile.values(:, 1, 1) = [
    420;
    600;
]-1e-3;
lp = {load_profile, load_profile};
lp{2}.values(:, 1, 1) = [590; 590];
ePg = { [ 370 500; [50 100]-1e-3; -lp{1}.values(:, 1, 1)' ], ...
        [ 500 500; [90 90]; -lp{2}.values(:, 1, 1)' ] };
eLMP = {[25 35; 25 35; 25 35], [30 30; 30 30; 30 30]};
eTLMP = {[25 35; 30 30; 25 35], [30 30; 30 30; 30 30]};
for w = 1:2
    profiles = getprofiles(lp{w});
    nt = size(profiles(1).values, 1);       % number of periods
    mdi = loadmd(mpc, nt, xgd, [], [], profiles);
    if w > 1
        mdi.InitialPg = ms.Pg(:, 1);
    end
    mdo = most(mdi, mpopt);
    ms = most_summary(mdo);

    t_is(mdo.results.success, 1, 12, sprintf('%s window %d : %s', t, w, 'success'));
    t_is(ms.Pg, ePg{w}, 12, sprintf('%s window %d : %s', t, w, 'Pg'));
    t_is(mdo.results.GenPrices, eLMP{w}, 5, sprintf('%s window %d : %s', t, w, 'LMP'));
    t_is(mdo.results.GenTLMP, eTLMP{w}, 5, sprintf('%s window %d : %s', t, w, 'TLMP'));

    if verbose
        fprintf('    ');
        for t = 1:nt
            fprintf('%-20s', sprintf('   T=%d', t));
        end
        fprintf('\n');
        for g = 1:2
            fprintf('G%d  ', g);
            for t = 1:nt
                fprintf('%-20s', ...
                    sprintf('(%.4g, %.4g, %.4g)', ms.Pg(g, t), ...
                                            mdo.results.GenPrices(g, t), ...
                                            mdo.results.GenTLMP(g, t) ));
            end
            fprintf('\n');
        end
    end
end


%%-----  TLMP Example 3  -----
t = 'example 3';
if verbose > 1
    fprintf('\n--------------------  TLMP Example 3  --------------------\n');
end
mpc = loadcase('t_case_tlmp_storage');
xgd_table.colnames = {
    'PositiveLoadFollowReservePrice', ...
            'PositiveLoadFollowReserveQuantity', ...
                'NegativeLoadFollowReservePrice', ...
                        'NegativeLoadFollowReserveQuantity', ...
};
xgd_table.data = [
    1e-6    40  1e-6    40;
    1e-6    60  1e-6    60;
    1e-6    Inf 1e-6    Inf;
];
load_profile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REP, ...
    'values', [] );
load_profile.values(:, 1, 1) = [
    199.9;
    219.1;
];% -1e-3;
lp = {load_profile, load_profile};
lp{2}.values(:, 1, 1) = [219.1; 118.8999];
%% StorageData
iess = 3;
mpc.iess = iess;
sv = 0;     %% cost/value of initial/residual stored energy
sd_table.OutEff     = 1;
sd_table.InEff      = 1;
sd_table.LossFactor = 0;
sd_table.rho        = 1;
sd_table.colnames = {
    'UnitIdx', ...
        'InitialStorage', ...
            'InitialStorageCost', ...
                'TerminalStoragePrice', ...
                    'MinStorageLevel', ...
                        'MaxStorageLevel', ...
                            'InitialStorageLowerBound', ...
                                'InitialStorageUpperBound', ...
};
sd_table.data = [
    iess    20  sv  sv  1   50  20  20;
];
ePg = { [ 199.9 200; -199.9 -219; 0 19 ], ...
        [ 159.9 119.9; -178.8999 -118.8999; 18.9999 -1.0001 ] };
eLMP = {[16.3 17.7; 16.3 17.7; 16.3 17.7], [15.7 13.6; 15.7 13.6; 15.7 13.6]};
eTLMP = {[16.3 17.7; 16.3 17.7; 16.3 17.7], [16.3 16.3; 17.7 11.6; 15.7 13.6]};
eTLMPs = {[14.3 15.7], [15.7 13.6]};
eSoC = {[20 1], [1.0001 2.0002]};
xgd = loadxgendata(xgd_table, mpc);
sd = loadstoragedata(sd_table, mpc.gen(iess, :));
for w = 1:2
    profiles = getprofiles(lp{w});
    nt = size(profiles(1).values, 1);       % number of periods
    mdi = loadmd(mpc, nt, xgd, sd, [], profiles);
    if w > 1
        mdi.InitialPg = ms.Pg(:, 1);
    end

    mdo = most(mdi, mpopt);
    ms = most_summary(mdo);

    t_is(mdo.results.success, 1, 12, sprintf('%s window %d : %s', t, w, 'success'));
    t_is(ms.Pg, ePg{w}, 5, sprintf('%s window %d : %s', t, w, 'Pg'));
    t_is(mdo.results.GenPrices, eLMP{w}, 5, sprintf('%s window %d : %s', t, w, 'LMP'));
    t_is(mdo.results.GenTLMP, eTLMP{w}, 5, sprintf('%s window %d : %s', t, w, 'GenTLMP'));
    t_is(mdo.results.StorageTLMPc, eTLMPs{w}, 5, sprintf('%s window %d : %s', t, w, 'StorageTLMPc'));
    t_is(mdo.results.StorageTLMPd, eTLMPs{w}, 5, sprintf('%s window %d : %s', t, w, 'StorageTLMPd'));
    t_is(ms.SoC, eSoC{w}, 5, sprintf('%s window %d : %s', t, w, 'SoC'));

    if verbose
        fprintf('    ');
        for t = 1:nt
            fprintf('%-20s', sprintf('   T=%d', t));
        end
        fprintf('\n');
        for g = 1:3
            fprintf('G%d  ', g);
            for t = 1:nt
                fprintf('%-20s', ...
                    sprintf('(%.4g, %.4g, %.4g)', ms.Pg(g, t), ...
                                            mdo.results.GenPrices(g, t), ...
                                            mdo.results.GenTLMP(g, t) ));
            end
            fprintf('\n');
        end
        for s = 1:ms.ns
            fprintf('S%d  ', s);
            for t = 1:nt
                fprintf('%-20s', ...
                    sprintf('(%.4g, %.4g, %.4g)', ms.SoC(s, t), ...
                                            mdo.results.StorageTLMPc(s, t), ...
                                            mdo.results.StorageTLMPd(s, t) ));
            end
            fprintf('\n');
        end
    end
end

%%-----  TLMP Example 4  -----
t = 'example 4';
% to match Cong's result, requires manually constraining storage to discharge
% at 0.1 in period 1
if verbose > 1
    fprintf('\n--------------------  TLMP Example 4  --------------------\n');
end
mpc = loadcase('t_case_tlmp_storage');
mpc.gencost(2, 5) = 1000;   %% make demand practically inelastic
xgd_table.colnames = {
    'PositiveLoadFollowReservePrice', ...
            'PositiveLoadFollowReserveQuantity', ...
                'NegativeLoadFollowReservePrice', ...
                        'NegativeLoadFollowReserveQuantity', ...
};
xgd_table.data = [
    1e-6    40  1e-6    40;
    1e-6    Inf 1e-6    Inf;
    1e-6    Inf 1e-6    Inf;
];
xgd_table.data = [
    1e-7    40  1e-7    40;
    1e-7    Inf 1e-7    Inf;
    1e-7    Inf 1e-7    Inf;
];
load_profile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REP, ...
    'values', [] );
load_profile.values(:, 1, 1) = [
    199.9;
    218.89999;
];% -1e-3;
lp = {load_profile, load_profile};
lp{2}.values(:, 1, 1) = [218.8999; 158.7];
%% StorageData
s_profile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TGEN, ...
    'rows', 3, ...
    'col', PMIN, ...
    'chgtype', CT_REP, ...
    'values', [] );
s_profile.values(:, 1, 1) = [
    0.1;
    -25;
];% -1e-3;
iess = 3;
mpc.iess = iess;
sv = 0;     %% cost/value of initial/residual stored energy
sd_table.OutEff     = 1;
sd_table.InEff      = 1;
sd_table.LossFactor = 0;
sd_table.rho        = 1;
sd_table.colnames = {
    'UnitIdx', ...
            'InitialStorage', ...
                'InitialStorageCost', ...
                    'TerminalStoragePrice', ...
                        'MinStorageLevel', ...
                            'MaxStorageLevel', ...
%                                 'InitialStorageLowerBound', ...
%                                     'InitialStorageUpperBound', ...
};
sd_table.data = [
    iess    20  sv  sv  1   50;
%     iess    20  sv  sv  1   50  20  20;
];
ePg = { [ 199.8 200; -199.9 -218.9; 0.1 18.9 ], ...
        [ 200 160; -218.9 -158.7; 18.9 -1.3 ] };
eLMP = {[16.3 16.3; 16.3 16.3; 16.3 16.3], [19 13.6; 19 13.6; 19 13.6]};
eTLMP = {[16.3 16.3; 16.3 16.3; 16.3 16.3], [16.3 16.3; 19 13.6; 19 13.6]};
eTLMPs = {[15.7 15.7], [15.7 13.6]};
eSoC = {[19.9 1], [1 2.3]};
xgd = loadxgendata(xgd_table, mpc);
sd = loadstoragedata(sd_table, mpc.gen(iess, :));
for w = 1:2
    profiles = getprofiles(lp{w});
    if w == 1   %% nudge storage dispatch to Cong's result
        profiles = getprofiles(profiles, s_profile);
    end
    nt = size(profiles(1).values, 1);       % number of periods
    mdi = loadmd(mpc, nt, xgd, sd, [], profiles);
    if w > 1
        mdi.InitialPg = ms.Pg(:, 1);
        mdi.Storage.InitialStorage = ms.SoC(:, 1);
    end

    mdo = most(mdi, mpopt);
    ms = most_summary(mdo);

    t_is(mdo.results.success, 1, 12, sprintf('%s window %d : %s', t, w, 'success'));
    t_is(ms.Pg, ePg{w}, 3.9, sprintf('%s window %d : %s', t, w, 'Pg'));
    t_is(mdo.results.GenPrices, eLMP{w}, 6, sprintf('%s window %d : %s', t, w, 'LMP'));
    t_is(mdo.results.GenTLMP, eTLMP{w}, 6, sprintf('%s window %d : %s', t, w, 'GenTLMP'));
    t_is(mdo.results.StorageTLMPc, eTLMPs{w}, 6, sprintf('%s window %d : %s', t, w, 'StorageTLMPc'));
    t_is(mdo.results.StorageTLMPd, eTLMPs{w}, 6, sprintf('%s window %d : %s', t, w, 'StorageTLMPd'));
    t_is(ms.SoC, eSoC{w}, 3.9, sprintf('%s window %d : %s', t, w, 'SoC'));

    if verbose
        fprintf('    ');
        for t = 1:nt
            fprintf('%-20s', sprintf('   T=%d', t));
        end
        fprintf('\n');
        for g = 1:3
            fprintf('G%d  ', g);
            for t = 1:nt
                fprintf('%-20s', ...
                    sprintf('(%.4g, %.4g, %.4g)', ms.Pg(g, t), ...
                                            mdo.results.GenPrices(g, t), ...
                                            mdo.results.GenTLMP(g, t) ));
            end
            fprintf('\n');
        end
        for s = 1:ms.ns
            fprintf('S%d  ', s);
            for t = 1:nt
                fprintf('%-20s', ...
                    sprintf('(%.4g, %.4g, %.4g)', ms.SoC(s, t), ...
                                            mdo.results.StorageTLMPc(s, t), ...
                                            mdo.results.StorageTLMPd(s, t) ));
            end
            fprintf('\n');
        end
    end
end

t_end

