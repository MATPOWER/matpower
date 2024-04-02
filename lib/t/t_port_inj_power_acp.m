function obj = t_port_inj_power_acp(quiet)
% t_port_inj_power_acp - Tests of mp.form_ac.port_inj_power derivatives wrt polar V.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

tc = struct( ...        %% test cases
        'name', {'1', '2'}, ...
        'ec', { @mp.nme_gizmo_acp, ...
                {   {@mp.nme_gen_acp_nln, 'mp.nme_gen'}, ...
                    {@mp.nme_load_acp_nln, 'mp.nme_load'}, ...
                    {@mp.nme_branch_acp_nln, 'mp.nme_branch'}, ...
                    {@mp.nme_shunt_acp_nln, 'mp.nme_shunt'}, ...
                    @mp.nme_gizmo_acp_nln }    } ...
    );

t_begin(87*length(tc), quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_gizmo';
mpopt = mpoption('out.all', 0, 'verbose', 0);
dmc = mp.dm_converter_mpc2().modify_element_classes(@mp.dmce_gizmo_mpc2).build();

for c = 1:length(tc)
    %% create network model object
    mpc = loadcase(casefile);
    dm = mp.data_model().modify_element_classes(@mp.dme_gizmo).build(mpc, dmc);
    ac = mp.net_model_acp().modify_element_classes(@mp.nme_gizmo_acp).build(dm);
    C = ac.C;
    D = ac.D;
    np = ac.np;
    nv = ac.nv/2;
    nz = ac.nz;
    A = [   C sparse(nv, nz);
            sparse(nz, np) D    ];
    A2 = [A sparse(nv+nz, np+nz); sparse(nv+nz, np+nz) A];

    %% other parameters
    dx = 1e-8;
    idx = randperm(np, fix(0.67*np))';
    lam = (1.5*rand(np, 1) + 0.5); k = randperm(np, fix(np/2)); lam(k) = -lam(k);
    e0 = zeros(np, 1);
    e1 = ones(np, 1);

    %% get ref bus index
    ref = find(mpc.bus(:, BUS_TYPE) == REF);

    %% construct initial system v1, v2, zr, zi, v_, z_, x_
    t = sprintf('%s : construct initial system v_, z_', tc(c).name);
    sv1 = ac.params_var('va');
    sv2 = ac.params_var('vm');
    szr = ac.params_var('zr');
    szi = ac.params_var('zi');

    %% randomize voltages a bit
    sv1 = sv1 + (0.6*rand(size(sv1)) - 0.3);    sv1(ref) = 0;
    sv2 = sv2 + (0.06*rand(size(sv1)) - 0.03);  sv2(ref) = 1;

    %% adjust values of z_
    szr(1) = 0.67;
    szi(1:3) = [0.1; 0.2; 0.3];

    %% initialize v_, z_, x_
    sv = sv2 .* exp(1j * sv1);
    sz = szr + 1j * szi;
    sx = [sv; sz];
    nx = length(sx);
    t_is(nx, nv+nz, 12, t);

    %%-----  tests using system voltages  -----
    t = sprintf('%s : ac.port_inj_power(x_) : ', tc(c).name);
    v10 = sv1; v20 = sv2; zr0 = szr; zi0 = szi; %% init w/ system v_, z_ components
    v0 = v20 .* exp(1j * v10);
    z0 = zr0 + 1j * zi0;
    x0 = [v0; z0];
    Nv = length(v0);
    Nz = length(z0);
    [S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0);   %% analytical

    %% check matrix input/output
    SS0 = ac.port_inj_power(x0*ones(1,Nv));
    t_is(SS0, S0*ones(1,Nv), 12, [t 'matrix input']);

    %% Sv1
    v_ = (v20*ones(1,Nv)) .* exp(1j * (v10*ones(1,Nv) + dx*eye(Nv,Nv)));
    z_ = (zr0 + 1j * zi0) * ones(1,Nv);
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_);
    num_Sv1b = (SS - SS0) / dx;
    t_is(full(Sv1), num_Sv1b, 6, [t 'Sv1']);

    %% Sv2
    v_ = (v20*ones(1,Nv) + dx*eye(Nv,Nv)) .* exp(1j * (v10*ones(1,Nv)));
    z_ = (zr0 + 1j * zi0) * ones(1,Nv);
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_);
    num_Sv2b = (SS - SS0) / dx;
    t_is(full(Sv2), num_Sv2b, 6, [t 'Sv2']);

    SS0 = ac.port_inj_power(x0*ones(1,Nz));

    %% Szr
    v_ = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
    z_ = (zr0 * ones(1,Nz) + dx*eye(Nz,Nz)) + 1j * (zi0 * ones(1,Nz));
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_);
    num_Szrb = (SS - SS0) / dx;
    t_is(full(Szr), num_Szrb, 6, [t 'Szr']);

    %% Szi
    v_ = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
    z_ = (zr0 * ones(1,Nz)) + 1j * (zi0 * ones(1,Nz) + dx*eye(Nz,Nz));
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_);
    num_Szib = (SS - SS0) / dx;
    t_is(full(Szi), num_Szib, 6, [t 'Szi']);

    t = sprintf('%s : ac.port_inj_power(x_, 1, idx) : ', tc(c).name);
    [iS0, iSv1, iSv2, iSzr, iSzi] = ac.port_inj_power(x0, 1, idx);
    t_is(iS0, S0(idx), 12, [t 'S0']);
    t_is(iSv1, Sv1(idx, :), 12, [t 'Sv1']);
    t_is(iSv2, Sv2(idx, :), 12, [t 'Sv2']);
    t_is(iSzr, Szr(idx, :), 12, [t 'Szr']);
    t_is(iSzi, Szi(idx, :), 12, [t 'Szi']);

    t = sprintf('%s : ac.port_inj_power_hess(x0, ek) == ac.p_i_p_h(x0, 1, 1, k) : ', tc(c).name);
    for k = 1:length(lam)
        ek = e0; ek(k) = 1;
        H1 = ac.port_inj_power_hess(x0, ek);
        H2 = ac.port_inj_power_hess(x0, 1, 1, k);
        t_is(H1, H2, 12, sprintf('%s%d', t, k));
    end

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam);
    HH = sparse(size(H, 1), size(H, 2));
    for k = 1:length(lam)
        ek = e0; ek(k) = 1;
        HH = HH + lam(k) * ac.port_inj_power_hess(x0, ek);
    end
    t_is(H, HH, 12, [t 'weighted sum indiv Hessians']);

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam, 1, idx) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam(idx), 1, idx);
    HH = sparse(size(H, 1), size(H, 2));
    for k = 1:length(idx)
        ek = e0; ek(idx(k)) = 1;
        HH = HH + lam(idx(k)) * ac.port_inj_power_hess(x0, ek);
    end
    t_is(H, HH, 12, [t 'weighted sum indiv Hessians']);

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam);
    [S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0);
    numH = zeros(2*nx, 2*nx);
    for k = 1:Nv
        v1 = v10; v1(k) = v1(k) + dx;
        v_ = v20 .* exp(1j * v1);
        x_ = [v_; z0];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_);
        numH(:, k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

        v2 = v20; v2(k) = v2(k) + dx;
        v_ = v2 .* exp(1j * v10);
        x_ = [v_; z0];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_);
        numH(:, Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
    end
    for k = 1:Nz
        z_ = zr0 + 1j * zi0;
        z_(k) = z_(k) + dx;
        x_ = [v0; z_];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_);
        numH(:, 2*Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

        z_ = zr0 + 1j * zi0;
        z_(k) = z_(k) + 1j*dx;
        x_ = [v0; z_];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_);
        numH(:, 2*Nv+Nz+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
    end
    t_is(full(H), numH, 5, [t 'numerical Hessian']);

    %%-----  tests using port voltages  -----
    t = sprintf('%s : ac.port_inj_power(x_, 0) : ', tc(c).name);
    v10 = C'*sv1; v20 = C'*sv2; zr0 = D'*szr; zi0 = D'*szi; %% init w/ port v_, z_ components
    v0 = v20 .* exp(1j * v10);
    z0 = zr0 + 1j * zi0;
    x0 = [v0; z0];
    Nv = length(v0);
    Nz = length(z0);
    [S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0, 0);    %% analytical

    %% check matrix input/output
    SS0 = ac.port_inj_power(x0*ones(1,Nv), 0);
    t_is(SS0, S0*ones(1,Nv), 12, [t 'matrix input']);

    %% Sv1
    v_ = (v20*ones(1,Nv)) .* exp(1j * (v10*ones(1,Nv) + dx*eye(Nv,Nv)));
    z_ = (zr0 + 1j * zi0) * ones(1,Nv);
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_, 0);
    num_Sv1b = (SS - SS0) / dx;
    t_is(full(Sv1), num_Sv1b, 6, [t 'Sv1']);

    %% Sv2
    v_ = (v20*ones(1,Nv) + dx*eye(Nv,Nv)) .* exp(1j * (v10*ones(1,Nv)));
    z_ = (zr0 + 1j * zi0) * ones(1,Nv);
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_, 0);
    num_Sv2b = (SS - SS0) / dx;
    t_is(full(Sv2), num_Sv2b, 6, [t 'Sv2']);

    SS0 = ac.port_inj_power(x0*ones(1,Nz), 0);

    %% Szr
    v_ = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
    z_ = (zr0 * ones(1,Nz) + dx*eye(Nz,Nz)) + 1j * (zi0 * ones(1,Nz));
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_, 0);
    num_Szrb = (SS - SS0) / dx;
    t_is(full(Szr), num_Szrb, 6, [t 'Szr']);

    %% Szi
    v_ = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
    z_ = (zr0 * ones(1,Nz)) + 1j * (zi0 * ones(1,Nz) + dx*eye(Nz,Nz));
    x_ = [v_; z_];
    SS = ac.port_inj_power(x_, 0);
    num_Szib = (SS - SS0) / dx;
    t_is(full(Szi), num_Szib, 6, [t 'Szi']);

    t = sprintf('%s : ac.port_inj_power(x_, 0, idx) : ', tc(c).name);
    [iS0, iSv1, iSv2, iSzr, iSzi] = ac.port_inj_power(x0, 0, idx);
    t_is(iS0, S0(idx), 12, [t 'S0']);
    t_is(iSv1, Sv1(idx, :), 12, [t 'Sv1']);
    t_is(iSv2, Sv2(idx, :), 12, [t 'Sv2']);
    t_is(iSzr, Szr(idx, :), 12, [t 'Szr']);
    t_is(iSzi, Szi(idx, :), 12, [t 'Szi']);

    t = sprintf('%s : ac.port_inj_power_hess(x0, ek, 0) == ac.p_i_p_h(x0, 1, 0, k) : ', tc(c).name);
    for k = 1:length(lam)
        ek = e0; ek(k) = 1;
        H1 = ac.port_inj_power_hess(x0, ek, 0);
        H2 = ac.port_inj_power_hess(x0, 1, 0, k);
        t_is(H1, H2, 12, sprintf('%s%d', t, k));
    end

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam, 0) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam, 0);
    HH = sparse(size(H, 1), size(H, 2));
    for k = 1:length(lam)
        ek = e0; ek(k) = 1;
        HH = HH + lam(k) * ac.port_inj_power_hess(x0, ek, 0);
    end
    t_is(H, HH, 12, [t 'weighted sum indiv Hessians']);

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam, 0, idx) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam(idx), 0, idx);
    HH = sparse(size(H, 1), size(H, 2));
    for k = 1:length(idx)
        ek = e0; ek(idx(k)) = 1;
        HH = HH + lam(idx(k)) * ac.port_inj_power_hess(x0, ek, 0);
    end
    t_is(H, HH, 12, [t 'weighted sum indiv Hessians']);

    t = sprintf('%s : ac.port_inj_power_hess(x_, lam, 0) : ', tc(c).name);
    H = ac.port_inj_power_hess(x0, lam, 0);
    [S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0, 0);
    Nx = 2*Nv+2*Nz;
    numH = zeros(Nx, Nx);
    for k = 1:Nv
        v1 = v10; v1(k) = v1(k) + dx;
        v_ = v20 .* exp(1j * v1);
        x_ = [v_; z0];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_, 0);
        numH(:, k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

        v2 = v20; v2(k) = v2(k) + dx;
        v_ = v2 .* exp(1j * v10);
        x_ = [v_; z0];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_, 0);
        numH(:, Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
    end
    for k = 1:Nz
        z_ = zr0 + 1j * zi0;
        z_(k) = z_(k) + dx;
        x_ = [v0; z_];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_, 0);
        numH(:, 2*Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

        z_ = zr0 + 1j * zi0;
        z_(k) = z_(k) + 1j*dx;
        x_ = [v0; z_];
        [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x_, 0);
        numH(:, 2*Nv+Nz+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
    end
    t_is(full(H), numH, 5, [t 'numerical Hessian']);
end

t_end;

if nargout
    obj = ac;
end
