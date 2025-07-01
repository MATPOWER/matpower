function mpc3p = convert_1p_to_3p(mpc, basekVA, basekV, freq, ishybrid)
    % convert_1p_to_3p - converts a single-phase MATPOWER case to an equivalent
    %   balanced three-phase case
    %   ::
    %
    %   MPC3P = CONVERT_1P_TO_3P(MPC, BASEKVA)
    %   MPC3P = CONVERT_1P_TO_3P(MPC, BASEKVA, BASEKV)
    %   MPC3P = CONVERT_1P_TO_3P(MPC, BASEKVA, BASEKV, FREQ)
    %   MPC3P = CONVERT_1P_TO_3P(MPC, BASEKVA, BASEKV, FREQ, ISHYBRID)
    %
    %   The conversion is performed taking into account the data format prototype
    %   for unbalanced three-phase elements included in MATPOWER 8.0., as
    %   well as subsequent additions regarding three-phase shunt elements and tap
    %   parameters for three-phase transformers.
    %
    %   This function is useful for creating new test cases for unbalanced
    %   networks. The core idea is to generate an equivalent balanced three-phase
    %   network, which users can then convert into an unbalanced one by manually
    %   adjusting loads, lines, transformers, among other elements.
    %
    %   Inputs (all optional except MPC):
    %       mpc : a MATPOWER case. It can be a struct with all the fields of
    %             a valid MATPOWER case or a string with the name of an mpc
    %       basekVA: a positive scalar denoting the base power in kVA for
    %                convertion to per-unit values of the three-phase system
    %       basekV : a vector denoting the base voltages in kV of all buses
    %                for convertion to per-unit values of the three-phase system
    %       freq : a positive scalar denoting the frequency in Hz of the system
    %       ishybrid : when set to 1, it returns the original mpc with
    %                  additional fields for the three-phase equivalent of
    %                  the network. The default is 0, which returns empty
    %                  fields for the single-phase network. The buslink
    %                  field is always empty.
    %
    %   Outputs:
    %       mpc3p (struct) : a struct with/without the fields of the single-phase
    %                        case depending on the values of ishybrid, plus
    %                        the additional fields for the equivalent
    %                        three-phase network. The resulting case ensures
    %                        a sequential numbering of the three-phase nodes
    %                        via the ext2int MATPOWER function.
    %
    %   If a power flow is run over the mpc3p provided by this function, the
    %   results are the three-phase balanced equivalent of the results of
    %   the single-phase original case.
    %
    %   Example:
    %
    %     run_pf('case10ba')
    %     mpc3p = convert_1p_to_3p('case10ba');
    %     run_pf(mpc3p, mpoption, 'mpx', mp.xt_3p)
    %
    % See also t_convert_1p_to_3p

    %   MATPOWER
    %   Copyright (c) 2019-2025, Power Systems Engineering Research Center (PSERC)
    %   by Wilson Gonzalez Vanegas, Universidad Nacional de Colombia Sede Manizales
    %
    %   This file is part of MATPOWER.
    %   Covered by the 3-clause BSD License (see LICENSE file for details).
    %   See https://matpower.org for more info.

    %% check inputs
    if nargin >= 1
        mpc = loadcase(mpc);
    else
        error('convert_1p_to_3p: A MATPOWER case must be provided \n')
    end

    BASE_KV = 10;            % Define index for base voltage in bus matrix
    change_base_kV  = 0;     % do not change base voltage by default
    change_base_kVA = 0;     % do not change base power by default
    if nargin < 5
        ishybrid = 0;
        if nargin < 4
            freq = 60;
            if nargin < 3
                basekV = mpc.bus(:, BASE_KV);
                if nargin < 2
                    basekVA = 1000 * mpc.baseMVA;
                    basekV_old = mpc.bus(:,BASE_KV);   % old base voltages
                else
                    change_base_kVA = 1;               % change base power
                    basekVA_old = 1000 * mpc.baseMVA;  % old base power
                    basekV_old = mpc.bus(:,BASE_KV);   % old base voltages
                end
            else
                change_base_kV = 1;                    % change base voltage
                change_base_kVA = 1;                   % change base power
                basekV_old = mpc.bus(:,BASE_KV);       % old base voltages
                basekVA_old = 1000 * mpc.baseMVA;      % old base power
            end
        end
    end

    if isstruct(mpc)
        flds = {'baseMVA', 'bus', 'branch', 'gen'};
        membs = ismember(fieldnames(mpc), flds);
        if sum(membs) ~= length(flds)
            error('convert_1p_to_3p: Not a valid MATPOWER case. Check for missing fields \n')
        end
    elseif isstring(mpc) || ischar(mpc)
        mpc = loadcase(mpc);
    else
        error('convert_1p_to_3p: Input ''mpc'' must be a struct or a string \n')
    end

    if ~isscalar(basekVA) || basekVA < 0
        error('convert_1p_to_3p: Input ''basekVA'' must be a positive scalar \n')
    end

    nb = size(mpc.bus, 1);
    if ~isvector(basekV) || numel(basekV) ~= nb || sum(basekV < 0) ~= 0
        error('convert_1p_to_3p: Input ''basekV'' must be a vector of dimension %d with positive voltage values in kV \n', nb)
    end

    if ~isscalar(freq) || freq < 0
        error('convert_1p_to_3p: Input ''freq'' must be a positive scalar, typically 50 or 60 \n')
    end

    %% define named indices into data matrices of single-phase case (mpc)
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
        VA, BASE_KV, ZONE, VMAX, VMIN] = idx_bus;
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
        ANGMIN, ANGMAX] = idx_brch;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

    %% classify branches into lines, transformers and general branches
    branch_is_xfmr = mpc.branch(:,TAP) ~= 0;
    line_has_series = ((mpc.branch(:,BR_R)) ~= 0 | (mpc.branch(:,BR_X)) ~= 0) & ~branch_is_xfmr;
    line_has_shunt = mpc.branch(:,BR_B) ~= 0 & ~branch_is_xfmr;
    branch_is_line = line_has_series | line_has_shunt;
    branch_is_general = mpc.branch(:,BR_B) ~= 0 & branch_is_xfmr;

    %% check for zero values in base voltages of three-phase nodes
    if ~all(basekV)
        error('convert_1p_to_3p: All nodal base voltages must be positive values in kV.')
    end

    %% check for presence of general branches
    if any(branch_is_general)
        error('convert_1p_to_3p: General branches with R or X, and B with an ideal tranformer are not supported. Please provide lines and transformers as individual elements in the branch matrix of the MATPOWER case.')
    end

    %% check for presence of phase shifters
    if any(mpc.branch(:, SHIFT))
        error('convert_1p_to_3p: Phase shifting transformer are not included yet in supported versions of the MATPOWER API for three-phase networks');
    end

    %% apply a reordering of bus indexing
    mpc = ext2int(mpc);     % facilitates work

    %% --- (1) bus3p: create data for three-phase buses
    nb = size(mpc.bus,1);                                   % number of single-phase buses
    bus3p = zeros(nb, 9);                                   % columns of bus3p are: busid (1), type (2), basekV (3), Vm1 (4), Vm2 (5), Vm3 (6), Va1 (7), Va2 (8), and Va3 (9)
    bus3p(:,1) = mpc.bus(:,BUS_I);                          % bus id
    bus3p(:,2) = mpc.bus(:,BUS_TYPE);                       % bus type
    if change_base_kV
        bus3p(mpc.bus(:,BUS_I),3) = basekV;                 % base voltage
        bus3p(mpc.bus(:,BUS_I),4) = ...                     % voltage magnitude for phase a [p.u]
            mpc.bus(:,VM) .* (basekV_old ./ basekV);
    else
        bus3p(mpc.bus(:,BUS_I),3) = mpc.bus(:,BASE_KV);     % base voltage
        bus3p(mpc.bus(:,BUS_I),4) = mpc.bus(:,VM);          % voltage magnitude for phase a [p.u]
    end
    bus3p(mpc.bus(:,BUS_I),5) = bus3p(mpc.bus(:,BUS_I),4);  % voltage magnitude for phase b [p.u]
    bus3p(mpc.bus(:,BUS_I),6) = bus3p(mpc.bus(:,BUS_I),4);  % voltage magnitude for phase c [p.u]
    bus3p(mpc.bus(:,BUS_I),7) = mpc.bus(:,VA);              % voltage angle for phase a [degrees]
    bus3p(mpc.bus(:,BUS_I),8) = mpc.bus(:,VA) - 120;        % voltage angle for phase b [degrees]
    bus3p(mpc.bus(:,BUS_I),9) = mpc.bus(:,VA) + 120;        % voltage angle for phase c [degrees]

    %% --- (2) xfmr3p: create data for three-phase transformers (Y-Y only, at the moment)
    if any(branch_is_xfmr)
        nx = sum(branch_is_xfmr);                           % number of three-phase transformers
        fbus_xfmr3p = mpc.branch(branch_is_xfmr, F_BUS);    % get buses at the from end of transformers
        tbus_xfmr3p = mpc.branch(branch_is_xfmr, T_BUS);    % get buses at the to end of transformers
        xfmr3p = ones(nx,9);                                % columns of xfmr3p are: xfid (1), fbus (2), tbus (3), status (4), R (5), X (6), basekVA (7), basekV (8), tap (9)
        xfmr3p(:,1) = (1:nx);                               % IDs of three-phase transformers
        xfmr3p(:,2) = fbus_xfmr3p;                          % from bus ID of three-phase transformers
        xfmr3p(:,3) = tbus_xfmr3p;                          % to bus ID of three-phase transformers
        r = mpc.branch(branch_is_xfmr, BR_R);
        x = mpc.branch(branch_is_xfmr, BR_X);
        if change_base_kV || change_base_kVA
            xfmr3p(:,5) = z_base_change(r, basekV_old(tbus_xfmr3p), basekVA_old, ...   % resistance of the three-phase transformers (p.u)
                                basekV(tbus_xfmr3p), basekVA);
            xfmr3p(:,6) = z_base_change(x, basekV_old(tbus_xfmr3p), basekVA_old, ...   % reactance of the three-phase transformers (p.u)
                                basekV(tbus_xfmr3p), basekVA);
        else
            xfmr3p(:,5) = r;                                                           % resistance of the three-phase transformers (p.u)
            xfmr3p(:,6) = x;                                                           % reactance of the three-phase transformers (p.u)
        end
        xfmr3p(:,7) = basekVA / 3;                          % base aparent power of the three-phase transformers
        xfmr3p(:,8) = basekV(fbus_xfmr3p) / sqrt(3);        % base voltage of the transformers
        xfmr3p(:,9) = mpc.branch(branch_is_xfmr, TAP);
    else
        xfmr3p = [];
    end

    %% --- (3) shunt3p: create data for three-phase shunts (wye-connected shunts, at the moment)
    bus_has_shunt = mpc.bus(:,GS) ~= 0 | mpc.bus(:,BS) ~= 0;
    if any(bus_has_shunt)
        ns = sum(bus_has_shunt);                            % numbers of shunts
        id_shunt = find(bus_has_shunt);                     % buses with shunts
        shunt3p = ones(ns,9);                               % columns of shunt3p are: shid (1), shbus (2), status (3), gs1 (4), gs2 (5), gs3 (6), bs1 (7), bs2 (8), bs3 (9)
        shunt3p(:,1) = (1:ns);                              % IDs of three-phase shunt
        shunt3p(:,2) = id_shunt;                            % bus id of three-phase shunt
        gs = mpc.bus(id_shunt, GS);
        bs = mpc.bus(id_shunt, BS);
        shunt3p(:,[4 5 6]) = repmat((1/3)*gs*1000, 1, 3);     % shunt conductances specified as nominal (@ vm = 1.0 p.u.) active power demand for all three phases [kW]
        shunt3p(:,[7 8 9]) = repmat((1/3)*bs*1000, 1, 3);     % shunt susceptances specified as nominal (@ vm = 1.0 p.u.) active power demand for all three phases [kVAr]
    else
        shunt3p = [];
    end

    %% --- (4) lc: create line constructions for building three-phase lines
    nl = sum(branch_is_line);                                      % number of three-phase lines
    fbus_branch = mpc.branch(:, F_BUS);
    tbus_branch = mpc.branch(:, T_BUS);

    % check for inconsistent base voltages at from/to ends of lines
    diff_ft_basekV = diff([basekV_old(fbus_branch) basekV_old(tbus_branch)],1,2);
    diff_ft_basekV = diff_ft_basekV~=0 & branch_is_line;
    if any(diff_ft_basekV)
        id_diff = find(diff_ft_basekV);
        error('convert_1p_to_3p: Base voltages at the from and to ends of a line must be equal. Check branches: %s',strjoin(cellstr(num2str(id_diff(:))), ', '));
    end

    fbus_line3p = mpc.branch(branch_is_line, F_BUS);                % get buses at from end of three-phase lines
    tbus_line3p = mpc.branch(branch_is_line, T_BUS);                % get buses at to end of three-phase lines
    baseY = mpc.baseMVA ./ (basekV_old(fbus_line3p)/sqrt(3)).^2;    % base admitance of lines using old bases (Ohm^-1)
    baseZ = 1./baseY;
    b = mpc.branch(branch_is_line, BR_B);
    r = mpc.branch(branch_is_line, BR_R);
    x = mpc.branch(branch_is_line, BR_X);
    R_lc = r .* baseZ * 3;                                          % resistance in Ohm for line construction
    X_lc = x .* baseZ * 3;                                          % reactance in Ohm for line construction
    C_lc = b .* baseY / 3 ./ (2*pi*freq*1e-9);                      % capacitance in nF for line construction
    line_construction = zeros(nl, 19);                              % columns of lc are: lcid (1), R11 (2), R21 (3), R31 (4), R22 (5), R32 (6), R33 (7), X11 (8), X21 (9), X31 (10), X22 (11), X32 (12), X33 (13), C11 (14), C21 (15), C31 (16), C22 (17), C32 (18), C33 (19)
    line_construction(:,1) = (1:nl);                                % line construction IDs (a lc for each line)
    line_construction(:,[2 5 7]) = repmat(R_lc, 1, 3);              % balanced diagonal matrix of resistances (Ohm/mile)
    line_construction(:,[8 11 13]) = repmat(X_lc, 1, 3);            % balanced diagonal matrix of reactances (Ohm/mile)
    line_construction(:,[14 17 19]) = repmat(C_lc, 1, 3);           % balanced diagonal matrix of capacitances (nF/mile)

    %% --- (5) line3p: create data for three-phase lines
    line3p = ones(nl, 6);                                   % columns of line3p are: brid (1), fbus (2), tbus (3), status (4), lcid (5), len (6)
    line3p(:,1) = (1:nl);                                   % IDs of three-phase lines
    line3p(:,2) = fbus_line3p;                              % from bus ID of three-phase lines
    line3p(:,3) = tbus_line3p;                              % to bus ID of three-phase lines
    line3p(:,5) = (1:nl);                                   % line construction IDs, each three-phase line has its own lc with length 1 mile

    %% --- (6) load3p: create data for three-phase loads (wye-connected loads, at the moment)
    id_Q_only = mpc.bus(:,PD) == 0 & mpc.bus(:,QD) ~= 0;                        % buses with reactive power-only loads
    id_load = (mpc.bus(:,PD) ~= 0 | mpc.bus(:,QD) ~= 0) & ~ id_Q_only;          % buses with active power loads
    nl_Q_only = sum(id_Q_only);
    nl = sum(id_load);                                                          % number of loads
    pd = mpc.bus(id_load,PD);
    qd = mpc.bus(id_load,QD);
    signPD = sign(pd); signPD(signPD==0) = 1;
    signQD = sign(qd); signQD(signQD==0) = 1;
    ldpf = signPD .* signQD .* cos(atan(qd ./ pd));                             % load power factor
    load3p = ones(nl+2*nl_Q_only, 9);                                           % columns of load3p are: ldid (1), ldbus (2), status (3), Pd1 (4) Pd2 (5), Pd3 (6), ldpf1 (7), ldpf2 (8), ldpf3 (9)
    load3p(:,1) = (1:size(load3p,1));                                           % IDs of thre-phase loads
    load3p(:,2) = [mpc.bus(id_load, BUS_I)                                      % bus id of three-phase loads
                   mpc.bus(id_Q_only, BUS_I)
                   mpc.bus(id_Q_only, BUS_I)];
    load3p(:,[4 5 6]) = [repmat((1/3)*mpc.bus(id_load, PD)*1000, 1, 3)          % active power demand for all three phases [kW]
                         repmat((1/3)*mpc.bus(id_Q_only, QD)*1000, 1, 3)        % dding a fictitious active power demand that is equal to the actual reactive power demand
                         -1*repmat((1/3)*mpc.bus(id_Q_only, QD)*1000, 1, 3)];   % adding the negative of the fictitious active power demand for exact modeling
    load3p(:,[7 8 9]) = [repmat(ldpf, 1, 3)                                     % load power factor for all three phases
                         repmat(sqrt(2)/2*ones(nl_Q_only,1), 1, 3)              % for fictitious active power demand, the angle is 45°, thus ldpf = sqrt(2)/2
                         -1*ones(nl_Q_only,3)];                                 % the power factor is -1 for the negative of the fictitious demand
    if nl_Q_only
        warning('convert_1p_to_3p: Reactive power-only loads were detected. To ensure exact modeling, additional negative active power loads were introduced at buses: %s.', strjoin(cellstr(num2str(mpc.bus(id_Q_only, BUS_I))), ', '))
    end

    %% --- (7) gen3p: create data for three-phase generators (wye-connected gens, at the moment)
    ng = size(mpc.gen, 1);                                              % number of generators
    gen3p = ones(ng, 12);                                               % columns of gen3p are: genid (1), gbus (2), status (3), Vg1 (4), Vg2 (5), Vg3 (6), Pg1 (7), Pg2 (8), Pg3 (9), Qg1 (10), Qg2 (11), Qg3 (12)
    gen3p(:,1) = (1:ng);                                                % IDs of three-phase generators
    gen3p(:,2) = mpc.gen(:,GEN_BUS);                                    % bus id of three-phase generators
    vg = mpc.gen(:,VG);
    if change_base_kV
        Vg = vg .* (basekV_old(gen3p(:,2)) ./ basekV(gen3p(:,2)));
        gen3p(:,[4 5 6]) = repmat(Vg, 1, 3);                            % voltage set point of three-phase generators
    else
        gen3p(:,[4 5 6]) = repmat(vg, 1, 3);                            % voltage set point of three-phase generators
    end
    gen3p(:,[7 8 9]) = repmat(1/3 * mpc.gen(:,PG) * 1000, 1 , 3);       % active power injection for all phases [kW]
    gen3p(:,[10 11 12]) = 0;                                            % reactive power injection for all phases [kVAr]

    %% create three-phase MATPOWER case (mpc3p)
    if ishybrid
        mpc3p = mpc;
    else
        mpc3p = struct('version', mpc.version, ...
                       'baseMVA', mpc.baseMVA, ...
                       'bus', [], ...
                       'gen', [], ...
                       'branch', [], ...
                       'gencost', []);
    end
    mpc3p.freq    = freq;
    mpc3p.basekVA = basekVA;
    mpc3p.bus3p   = bus3p;
    mpc3p.line3p  = line3p;
    mpc3p.xfmr3p  = xfmr3p;
    mpc3p.shunt3p = shunt3p;
    mpc3p.load3p  = load3p;
    mpc3p.gen3p   = gen3p;
    mpc3p.lc      = line_construction;
    mpc3p.buslink = [];
end

function z_new = z_base_change(z_old, basekV_old, basekVA_old, basekV_new, basekVA_new)
    z_new = z_old * (basekVA_new / basekVA_old) .* (basekV_old.^2) ./ (basekV_new.^2);
end