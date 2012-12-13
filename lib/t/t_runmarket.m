function t_runmarket(quiet)
%T_RUNMARKET  Tests for code in RUNMKT, SMARTMKT AND AUCTION.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2005-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

if nargin < 1
    quiet = 0;
end

n_tests = 20;

t_begin(n_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

if ~have_fcn('smartmarket')
    t_skip(n_tests, 'smartmarket code not available');
else
    mpc = loadcase('t_auction_case');
    
    mpopt = mpoption('OPF_ALG', 560, 'OUT_ALL_LIM', 1, 'OUT_BRANCH', 0, 'OUT_SYS_SUM', 0, 'OUT_ALL', 0, 'VERBOSE', 1);
    % mpopt = mpoption('OUT_GEN', 1, 'OUT_BRANCH', 0, 'OUT_SYS_SUM', 0);
    
    offers.P.qty = [
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
    ];
    offers.P.prc = [
        20 50 60;
        20 40 70;
        20 42 80;
        20 44 90;
        20 46 75;
        20 48 60;
    ];
    bids.P.qty = [
        10 10 10;
        10 10 10;
        10 10 10;
    ];
    bids.P.prc = [
        100 70 60;
%         100 64.3 20;
%         100 30.64545 0;
        100 50 20;
        100 60 50;
    ];
    
    offers.Q.qty = [ 60; 60; 60; 60; 60; 60; 0; 0; 0 ];
    offers.Q.prc = [ 0; 0; 0; 0; 0; 3; 0; 0; 0 ];
    bids.Q.qty = [ 15; 15; 15; 15; 15; 15; 15; 12; 7.5 ];
%     bids.Q.prc = [ 0; 0; 0; 0; 0; 0; 0; 83.9056; 0 ];
    bids.Q.prc = [ 0; 0; 0; 0; 0; 0; 0; 20; 0 ];
    
    t = 'marginal Q offer, marginal PQ bid, auction_type = 5';
    mkt = struct(   'auction_type', 5, ...
                    't' , [], ...
                    'u0', [], ...
                    'lim', []   );
    [r, co, cb, f, dispatch, success, et] = runmarket(mpc, offers, bids, mkt, mpopt);
    co5 = co;
    cb5 = cb;
    
%     [ co.P.qty co.P.prc ]
%     [ cb.P.qty cb.P.prc ]
%     [ co.Q.qty co.Q.prc ]
%     [ cb.Q.qty cb.Q.prc ]
    
    i2e = r.bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(r.bus, 1))';
    G = find( ~isload(r.gen) );   %% real generators
    L = find(  isload(r.gen) );   %% dispatchable loads
    Gbus = e2i(r.gen(G,GEN_BUS));
    Lbus = e2i(r.gen(L,GEN_BUS));
    
    t_is( co.P.qty, ones(6, 1) * [12 24 0], 2, [t ' : gen P quantities'] );
    t_is( co.P.prc(1,:), 50.1578*ones(1,3), 3, [t ' : gen 1 P prices'] );
    t_is( cb.P.qty, [10 10 10; 10 0.196 0; 10 10 0], 2, [t ' : load P quantities'] );
    t_is( cb.P.prc(2,:), 56.9853*ones(1,3), 4, [t ' : load 2 P price'] );
    t_is( co.P.prc(:,1), r.bus(Gbus, LAM_P), 8, [t ' : gen P prices'] );
    t_is( cb.P.prc(:,1), r.bus(Lbus, LAM_P), 8, [t ' : load P prices'] );
    
    t_is( co.Q.qty, [4.2722; 11.3723; 14.1472; 22.8939; 36.7886; 12.3375; 0; 0; 0], 2, [t ' : Q offer quantities'] );
    t_is( co.Q.prc, [0;0;0;0;0;3; 0.4861; 2.5367; 1.3763], 4, [t ' : Q offer prices'] );
    t_is( cb.Q.qty, [0;0;0;0;0;0; 15; 4.0785; 5], 2, [t ' : Q bid quantities'] );
    t_is( cb.Q.prc, [0;0;0;0;0;3; 0.4861; 2.5367; 1.3763], 4, [t ' : Q bid prices'] );
    t_is( co.Q.prc, r.bus([Gbus; Lbus], LAM_Q), 8, [t ' : Q offer prices'] );
    t_is( cb.Q.prc, co.Q.prc, 8, [t ' : Q bid prices'] );
    
    t = 'marginal Q offer, marginal PQ bid, auction_type = 0';
    mkt.auction_type = 0;
    [r, co, cb, f, dispatch, success, et] = runmarket(mpc, offers, bids, mkt, mpopt);
    t_is( co.P.qty, co5.P.qty, 8, [t ' : gen P quantities'] );
    t_is( cb.P.qty, cb5.P.qty, 8, [t ' : load P quantities'] );
    t_is( co.P.prc, offers.P.prc, 8, [t ' : gen P prices'] );
    t_is( cb.P.prc, bids.P.prc, 8, [t ' : load P prices'] );
    
    t_is( co.Q.qty, co5.Q.qty, 8, [t ' : gen Q quantities'] );
    t_is( cb.Q.qty, cb5.Q.qty, 8, [t ' : load Q quantities'] );
    t_is( co.Q.prc, offers.Q.prc, 8, [t ' : gen Q prices'] );
    t_is( cb.Q.prc, bids.Q.prc, 8, [t ' : load Q prices'] );
end

t_end;
