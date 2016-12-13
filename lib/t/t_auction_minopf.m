function t_auction_minopf(quiet)
%T_AUCTION_MINOPF  Tests for code in auction.m, using MINOPF solver.

%   MATPOWER
%   Copyright (c) 2004-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin < 1
    quiet = 0;
end

n_tests = 183;

t_begin(n_tests, quiet);

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

if ~have_fcn('smartmarket')
    t_skip(n_tests, 'smartmarket code not available');
elseif ~have_fcn('minopf')
    t_skip(n_tests, 't_auction_minopf requires MINOPF');
else
    mpopt = mpoption('opf.ac.solver', 'MINOPF', 'out.lim.all', 1, 'out.branch', 0, 'out.sys_sum', 0, 'out.all', 0, 'verbose', 0);
    q = [
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
        12 24 24; 
        10 10 10;
        10 10 10;
        10 10 10;
    ];

    %%-----  one offer block marginal @ $50  -----
    p = [
        20 50 60;
        20 40 70;
        20 42 80;
        20 44 90;
        20 46 75;
        20 48 60;
        100 70 60;
        100 50 20;
        100 60 50;
    ];

    t = 'one marginal offer @ $50, auction_type = 5';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1150, 100, [], [], mpopt);
    cq5 = cq;
    cp5 = cp;
    i2e = bus(:, BUS_I);
    e2i = sparse(max(i2e), 1);
    e2i(i2e) = (1:size(bus, 1))';
    G = find( ~isload(gen) );   %% real generators
    L = find(  isload(gen) );   %% dispatchable loads
    Gbus = e2i(gen(G,GEN_BUS));
    Lbus = e2i(gen(L,GEN_BUS));
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));

    t_is( cq(G(1),2:3), [23.32 0], 2, t );
    t_is( cp(G(1),:), 50, 4, t );
    t_is( cq(L(2),1:2), [10 0], 2, t );
    t_is( cp(L(2),:), 54.0312, 4, t );
    t_is( cp(G,1), bus(Gbus, LAM_P), 8, [t ' : gen prices'] );
    t_is( cp(L,1), bus(Lbus, LAM_P) + Qfudge(L,1), 8, [t ' : load prices'] );

    lao_X = p(G(1),2)/bus(Gbus(1), LAM_P);
    fro_X = p(G(6),3)/bus(Gbus(6), LAM_P);
    lab_X = p(L(3),2)/(bus(Lbus(3), LAM_P) + Qfudge(L(3),1));
    frb_X = p(L(2),2)/(bus(Lbus(2), LAM_P) + Qfudge(L(2),1));

    t_is( lao_X, 1, 4, 'lao_X');
    t_is( fro_X, 1.1324, 4, 'fro_X');
    t_is( lab_X, 1.0787, 4, 'lab_X');
    t_is( frb_X, 0.9254, 4, 'frb_X');

    t = 'one marginal offer @ $50, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1110, 100, [], [], mpopt);
    cp1 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 8, [t ' : prices'] );

    t = 'one marginal offer @ $50, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1120, 100, [], [], mpopt);
    cp2 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:)*fro_X, 8, [t ' : gen prices'] );
    t_is( cp(L(1:2),:), cp5(L(1:2),:)*fro_X, 8, [t ' : load 1,2 prices'] );
    t_is( cp(L(3),:), 60, 5, [t ' : load 3 price'] );   %% clipped by accepted bid

    t = 'one marginal offer @ $50, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1130, 100, [], [], mpopt);
    cp3 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5*lab_X, 8, [t ' : prices'] );

    t = 'one marginal offer @ $50, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1140, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), p(G(1),2), 8, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp5(G(2:6),:)*frb_X, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp5(L,:)*frb_X, 8, [t ' : load prices'] );

    t = 'one marginal offer @ $50, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp3, 8, [t ' : prices'] );
    p2 = p;
    p2(L,:) = [ 100 100 100;
                100   0   0;
                100 100   0 ];
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 5, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:)*fro_X, 4, [t ' : gen prices'] );
    t_is( cp(L,:), cp5(L,:)*fro_X, 4, [t ' : load prices'] ); %% load 3 not clipped as in FRO

    t = 'one marginal offer @ $50, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1170, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5 * (lao_X+lab_X)/2, 8, [t ' : prices'] );
    t_is( cp, (cp1 + cp3) / 2, 8, [t ' : prices'] );

    t = 'one marginal offer @ $50, auction_type = 8';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1180, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp1(G,:), 8, [t ' : gen prices'] );
    t_is( cp(L,:), cp3(L,:), 8, [t ' : load prices'] );

    t = 'one marginal offer @ $50, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1100, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, p, 8, [t ' : prices'] );


    %%-----  one bid block marginal @ $55  -----
    p(L(2),2) = 55;
    t = 'one marginal bid @ $55, auction_type = 5';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1150, 100, [], [], mpopt);
    cq5 = cq;
    cp5 = cp;
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));

    t_is( cq(G(1),2:3), [24 0], 2, t );
    t_is( cp(G(1),:), 50.016, 3, t );
    t_is( cq(L(2),1:2), [10 0.63], 2, t );
    t_is( cp(L(2),:), 55, 4, t );
    t_is( cp(G,1), bus(Gbus, LAM_P), 8, [t ' : gen prices'] );
    t_is( cp(L,1), bus(Lbus, LAM_P) + Qfudge(L,1), 8, [t ' : load prices'] );

    lao_X = p(G(1),2)/bus(Gbus(1), LAM_P);
    fro_X = p(G(6),3)/bus(Gbus(6), LAM_P);
    lab_X = p(L(2),2)/(bus(Lbus(2), LAM_P) + Qfudge(L(2),1));
    frb_X = p(L(3),3)/(bus(Lbus(3), LAM_P) + Qfudge(L(3),1));

    t_is( lao_X, 0.9997, 4, 'lao_X');
    t_is( fro_X, 1.1111, 4, 'fro_X');
    t_is( lab_X, 1, 4, 'lab_X');
    t_is( frb_X, 0.8960, 4, 'frb_X');

    t = 'one marginal bid @ $55, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1110, 100, [], [], mpopt);
    cp1 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5*lao_X, 8, [t ' : prices'] );

    t = 'one marginal bid @ $55, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1120, 100, [], [], mpopt);
    cp2 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:)*fro_X, 8, [t ' : gen prices'] );
    t_is( cp(L(1),:), cp5(L(1),:)*fro_X, 8, [t ' : load 1 price'] );
    t_is( cp(L(2),:), 55, 5, [t ' : load 2 price'] );
    t_is( cp(L(3),:), 60, 5, [t ' : load 3 price'] );

    t = 'one marginal bid @ $55, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1130, 100, [], [], mpopt);
    cp3 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 8, [t ' : prices'] );

    t = 'one marginal bid @ $55, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1140, 100, [], [], mpopt);
    cp4 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 50, 5, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp5(G(2:6),:)*frb_X, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp5(L,:)*frb_X, 8, [t ' : load prices'] );

    t = 'one marginal bid @ $55, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp1, 8, [t ' : prices'] );

    p2 = p;
    p2(G,:) = [ 0 0 100;
                0 0 100;
                0 0 100;
                0 0 100;
                0 0 100;
                0 0 100 ];
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:)*frb_X, 4, [t ' : gen prices'] );  %% gen 1, not clipped this time
    t_is( cp(L,:), cp4(L,:), 4, [t ' : load prices'] );

    t = 'one marginal bid @ $55, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1170, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5 * (lao_X+lab_X)/2, 8, [t ' : prices'] );
    t_is( cp, (cp1 + cp3) / 2, 8, [t ' : prices'] );

    t = 'one marginal bid @ $55, auction_type = 8';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1180, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp1(G,:), 8, [t ' : gen prices'] );
    t_is( cp(L,:), cp3(L,:), 8, [t ' : load prices'] );

    t = 'one marginal bid @ $55, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1100, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, p, 8, [t ' : prices'] );


    %%-----  one bid block marginal @ $54.50 and one offer block marginal @ $50  -----
    p(L(2),2) = 54.5;
    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 5';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1150, 100, [], [], mpopt);
    cq5 = cq;
    cp5 = cp;
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));

    t_is( cq(G(1),2:3), [23.74 0], 2, t );
    t_is( cp(G(1),:), 50, 4, t );
    t_is( cq(L(2),1:2), [10 0.39], 2, t );
    t_is( cp(L(2),:), 54.5, 4, t );
    t_is( cp(G,1), bus(Gbus, LAM_P), 8, [t ' : gen prices'] );
    t_is( cp(L,1), bus(Lbus, LAM_P) + Qfudge(L,1), 8, [t ' : load prices'] );

    lao_X = p(G(1),2)/bus(Gbus(1), LAM_P);
    fro_X = p(G(6),3)/bus(Gbus(6), LAM_P);
    lab_X = p(L(2),2)/(bus(Lbus(2), LAM_P) + Qfudge(L(2),1));
    frb_X = p(L(3),3)/(bus(Lbus(3), LAM_P) + Qfudge(L(3),1));

    t_is( lao_X, 1, 4, 'lao_X');
    t_is( fro_X, 1.1221, 4, 'fro_X');
    t_is( lab_X, 1, 4, 'lab_X');
    t_is( frb_X, 0.8976, 4, 'frb_X');

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1110, 100, [], [], mpopt);
    cp1 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 4, [t ' : prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1120, 100, [], [], mpopt);
    cp2 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:)*fro_X, 5, [t ' : gen prices'] );
    t_is( cp(L(1),:), cp5(L(1),:)*fro_X, 5, [t ' : load 1 price'] );
    t_is( cp(L(2),:), 54.5, 5, [t ' : load 2 price'] );
    t_is( cp(L(3),:), 60, 5, [t ' : load 3 price'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1130, 100, [], [], mpopt);
    cp3 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 8, [t ' : prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1140, 100, [], [], mpopt);
    cp4 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 50, 5, [t ' : gen 1 price'] );
    t_is( cp(G(2:5),:), cp5(G(2:5),:)*frb_X, 8, [t ' : gen 2-5 prices'] );
    t_is( cp(G(6),:), 48, 5, [t ' : gen 6 price'] );
    t_is( cp(L,:), cp5(L,:)*frb_X, 8, [t ' : load prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 4, [t ' : prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1170, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 4, [t ' : prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 8';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1180, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 4, [t ' : prices'] );

    t = 'marginal offer @ $50, bid @ $54.50, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1100, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, p, 8, [t ' : prices'] );


    %%-----  gen 1 at Pmin, load 3 block 2 marginal @ $60  -----
    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 5';
    p(L(2),2) = 50;     %% undo previous change
    p2 = p;
    p2(G(1),2:3) = [65 65];
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1150, 100, [], [], mpopt);
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));

    t_is( cp(G(1),:), 65, 4, [t ' : gen 1 price'] );
    t_is( cp(G(2),:), 54.2976, 4, [t ' : gen 2 price'] );
    cq5 = cq;
    cp5 = cp;
    cp_lam = cp5;
    cp_lam(1,:) = bus(Gbus(1), LAM_P);  %% unclipped

    lao_X = p2(G(6),2)/bus(Gbus(6), LAM_P);
    fro_X = p2(G(6),3)/bus(Gbus(6), LAM_P);
    lab_X = p2(L(3),2)/(bus(Lbus(3), LAM_P) + Qfudge(L(3),1));
    frb_X = p2(L(2),2)/(bus(Lbus(2), LAM_P) + Qfudge(L(2),1));

    t_is( lao_X, 0.8389, 4, 'lao_X');
    t_is( fro_X, 1.0487, 4, 'fro_X');
    t_is( lab_X, 1, 4, 'lab_X');
    t_is( frb_X, 0.8569, 4, 'frb_X');

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1110, 100, [], [], mpopt);
    cp1 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 65, 8, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp_lam(G(2:6),:)*lao_X, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp_lam(L,:)*lao_X, 8, [t ' : load prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1120, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 65, 8, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp_lam(G(2:6),:)*fro_X, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L(1:2),:), cp_lam(L(1:2),:)*fro_X, 8, [t ' : load 1-2 prices'] );
    t_is( cp(L(3),:), 60, 8, [t ' : load 3 price'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1130, 100, [], [], mpopt);
    cp3 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 65, 8, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp_lam(G(2:6),:), 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp_lam(L,:), 8, [t ' : load prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1140, 100, [], [], mpopt);
    cp4 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 65, 5, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp5(G(2:6),:)*frb_X, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp5(L,:)*frb_X, 8, [t ' : load prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp4, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1170, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G(1),:), 65, 4, [t ' : gen 1 price'] );
    t_is( cp(G(2:6),:), cp_lam(G(2:6),:) * (lao_X+lab_X)/2, 8, [t ' : gen 2-6 prices'] );
    t_is( cp(L,:), cp_lam(L,:) * (lao_X+lab_X)/2, 8, [t ' : load prices'] );
    t_is( cp, (cp1 + cp3) / 2, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 8';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1180, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp1(G,:), 8, [t ' : prices'] );
    t_is( cp(L,:), cp3(L,:), 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal bid @ $60, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1100, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, p2, 8, [t ' : prices'] );


    %%-----  gen 1 at Pmin, gen 6 block 3 marginal @ $60  -----
    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 5';
    p2(L,:) = [ 100 100 100;
                100   0   0;
                100 100   0 ];
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1150, 100, [], [], mpopt);
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));

    t_is( cp(G(1),:), 65, 4, [t ' : gen 1 price'] );
    t_is( cp(G(2),:), 57.1616, 4, [t ' : gen 2 price'] );
    cq5 = cq;
    cp5 = cp;
    cp_lam = cp5;
    cp_lam(1,:) = bus(Gbus(1), LAM_P);  %% unclipped

    lao_X = p2(G(6),3)/bus(Gbus(6), LAM_P);
    fro_X = p2(G(1),3)/bus(Gbus(1), LAM_P);
    lab_X = p2(L(3),2)/(bus(Lbus(3), LAM_P) + Qfudge(L(3),1));
    frb_X = p2(L(2),2)/(bus(Lbus(2), LAM_P) + Qfudge(L(2),1));

    t_is( lao_X, 1, 4, 'lao_X');
    t_is( fro_X, 1.1425, 4, 'fro_X');
    t_is( lab_X, 1.5813, 4, 'lab_X');
    t_is( frb_X, 0, 4, 'frb_X');

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1110, 100, [], [], mpopt);
    cp1 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp5, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1120, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp_lam*fro_X, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1130, 100, [], [], mpopt);
    cp3 = cp;
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp_lam*lab_X, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1140, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,1), [65;40;42;44;46;60], 4, [t ' : gen prices'] );
    t_is( cp(L,:), cp_lam(L,:)*frb_X, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1160, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp_lam*fro_X, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1170, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, cp_lam * (lao_X+lab_X)/2, 8, [t ' : prices'] );
    t_is( cp, (cp_lam + cp3) / 2, 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 8';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1180, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp(G,:), cp5(G,:), 8, [t ' : prices'] );
    t_is( cp(L,:), cp3(L,:), 8, [t ' : prices'] );

    t = 'gen 1 @ Pmin, marginal offer @ $60, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p2, 1100, 100, [], [], mpopt);
    t_is( cq, cq5, 8, [t ' : quantities'] );
    t_is( cp, p2, 8, [t ' : prices'] );


    %%-----  gen 2 decommitted, one offer block marginal @ $60  -----
    p(G(2),:) = p(G(2),:) + 100;

    t = 'price of decommited gen, auction_type = 5';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1150, 200, [], [], mpopt);
    cp5 = cp;
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));
    t_is(sum(cq(2,:)), 0, 8, t);
    t_is(cp(2,1), 59.194, 3, t);

%     Xo = p(1:6, :) ./ (diag(bus(Gbus, LAM_P)) * ones(size(p(G,:))));
%     ao = (cq(1:6, :) ~= 0);
%     ro = (cq(1:6, :) == 0);
%     Xb = p(7:9, :) ./ (diag(bus(Lbus, LAM_P) + gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:))));
%     ab = (cq(7:9, :) ~= 0);
%     rb = (cq(7:9, :) == 0);
%     aXo = ao .* Xo
%     rXo = ro .* Xo
%     aXb = ab .* Xb
%     rXb = rb .* Xb

    lao_X = p(G(6),3)/bus(Gbus(6), LAM_P);
    fro_X = p(G(1),3)/bus(Gbus(1), LAM_P);
    lab_X = p(L(1),2)/(bus(Lbus(1), LAM_P) + Qfudge(L(1),1));
    frb_X = p(L(1),3)/(bus(Lbus(1), LAM_P) + Qfudge(L(1),1));

    t_is( lao_X, 1, 4, 'lao_X');
    t_is( fro_X, 1.0212, 4, 'fro_X');
    t_is( lab_X, 1.1649, 4, 'lab_X');
    t_is( frb_X, 0.9985, 4, 'frb_X');

    t = 'price of decommited gen, auction_type = 1';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1110, 200, [], [], mpopt);
    t_is(cp(2,1), 59.194, 3, t);

    t = 'price of decommited gen, auction_type = 2';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1120, 200, [], [], mpopt);
    t_is(cp(2,1), cp5(2,1)*fro_X, 3, t);

    t = 'price of decommited gen, auction_type = 3';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1130, 200, [], [], mpopt);
    t_is(cp(2,1), cp5(2,1)*lab_X, 3, t);

    t = 'price of decommited gen, auction_type = 4';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1140, 200, [], [], mpopt);
    t_is(cp(2,1), cp5(2,1)*frb_X, 3, t);

    t = 'price of decommited gen, auction_type = 6';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1160, 200, [], [], mpopt);
    t_is(cp(2,1), cp5(2,1)*fro_X, 3, t);

    t = 'price of decommited gen, auction_type = 7';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1170, 200, [], [], mpopt);
    t_is(cp(2,1), cp5(2,1)*(lao_X+lab_X)/2, 3, t);

    t = 'price of decommited gen, auction_type = 0';
    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1100, 200, [], [], mpopt);
    t_is(cp(2,1), 120, 3, t);

    t = 'single block, marginal offer @ $50, auction_type = 5';
    q = [
        60; 
        36; 
        36; 
        36; 
        36; 
        36; 
        30;
        10;
        20;
    ];

    p = [
        50;
        40;
        42;
        44;
        46;
        48;
        100;
        100;
        100;
    ];

    [MVAbase, cq, cp, bus, gen, gencost, branch, f, dispatch, success, et] = ...
        runmkt('t_auction_case', q, p, 1150, 100, [], [], mpopt);
    t_is( cq(G(1)), 35.32, 2, t );
    t_is( cq(G(2:6)), q(G(2:6)), 8, [t ' : gen qtys'] ); 
    t_is( cp(G(1)), 50, 4, t );
    t_is( cq(L), q(L), 8, [t ' : load qtys'] );
    t_is( cp(L(2),:), 54.03, 2, t );
    t_is( cp(G), bus(Gbus, LAM_P), 8, [t ' : gen prices'] );
    Qfudge =  zeros(size(p));
    Qfudge(L,:) = diag(gen(L,QG) ./ gen(L,PG) .* bus(Lbus, LAM_Q)) * ones(size(p(L,:)));
    t_is( cp(L), bus(Lbus, LAM_P) + Qfudge(L,1), 8, [t ' : load prices'] );
end

t_end;
