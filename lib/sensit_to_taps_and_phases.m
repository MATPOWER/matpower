function [ sensitivity_of_H_to_U_taps_and_phasess ] = sensit_to_taps_and_phases( results, tap_changers, phase_shifters, Jacobian , tip_vozlisc)
%sensit_to_taps_and_phases  Calculates the first orded sensitivities of 
%regulated active powers and nodal voltages to taps and phases.
%   [sensitivity_of_H_to_U_taps_and_phasess] = sensit_to_taps_and_phases(results, tap_changers, phase_shifters, Jacobian , node_types)
%
%   Inputs :
%       results: results from a runpf(), interanly idnexed (use ext2int())
%       tap_changers: tap changers data as in taps_and_phases_analysis()
%       phase_shifters: phase shifter data as in taps_and_phases_analysis()
%       Jacobian: power missmatch Jacobian matrix with the appropriate
%       ordering, rows go P1, Q1, P2, Q2...., columns go V1, delta1, V2,
%       delta2..... 
%       tip_vozlisc: node types, 2 for PQ, 1 for PV, same as nubmer of
%       equations for each node, used for indexing
%
%   Output :
%       sensitivity_of_H_to_U_taps_and_phasess : first order sensitivities 
%       of p.u. changes in nodal votlages and line active powers to changes
%       in calculated taps and phases of transformers.
%
%   by Gorazd Bone, Faculty of Electrical Engineering, Ljubljana 

ind_node_eq = tip_vozlisc * 0;
for k=1:size(tip_vozlisc,1)
    ind_node_eq(k) = sum(tip_vozlisc(1:k-1)) + 1;
end

stevilo_tap_ch = size(tap_changers,1);
stevilo_pha_sh = size(phase_shifters,1);

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;  %#ok<*NASGU,*ASGLU>
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% Define Partial derivative matrices
% G is the LF problem power injection missmatch form (sum of nodal powers)
GpoU = zeros(sum(tip_vozlisc) , stevilo_tap_ch + stevilo_pha_sh); % derivative of LF problem w.r.t. taps & phases
GpoX = Jacobian; % derivative of LF problem w.r.t. voltages and angles
HpoU = zeros(stevilo_tap_ch + stevilo_pha_sh); %derivative of criterion function w.r.t. taps & phases
HpoX = zeros(stevilo_tap_ch + stevilo_pha_sh , sum(tip_vozlisc));%derivative of criterion function w.r.t. voltages and angles

for k=1:stevilo_tap_ch % partial derivatives for taps
    veja = tap_changers(k,1);%tap branch
    reg_voz = tap_changers(k,2); %ragulated node
    voz1 = results.branch(veja,F_BUS);
    voz2 = results.branch(veja,T_BUS);
    indeks_P1 = ind_node_eq(voz1); %active power and volt. amp index for f node
    indeks_P2 = ind_node_eq(voz2); %active power and volt. amp index for t node
    indeks_reg_voz = ind_node_eq(reg_voz);%volt. amp index for reg. node
    
    Uvozi = results.bus(voz1,VM);
    Uvozj = results.bus(voz2,VM);
    di = results.bus(voz1,VA)/180*pi;
    dj = results.bus(voz2,VA)/180*pi;
    Z = results.branch(veja,BR_R) + 1i*results.branch(veja,BR_X);
    Bsh = results.branch(veja,BR_B);
    G = real(1/Z);
    B = imag(1/Z);
    tap = results.branch(veja,TAP);
    ph = results.branch(veja,SHIFT)/180*pi;
    
    if tip_vozlisc(reg_voz) == 1 % tap changer regulating PV node
        error_string =  strcat(2 , 'Tap changer' , int2str(k) , 'regulating the voltage of a PV node');
        errordlg(error_string );
        return
    end
    
    dPijdtap = (Uvozi*(-2*G*Uvozi+G*tap*Uvozj*cos(di-dj-ph)+B*tap*Uvozj*sin(di-dj-ph)))/tap^3;
    dPjidtap = (Uvozi*Uvozj*(G*cos(di-dj-ph)-B*sin(di-dj-ph)))/tap^2;
    dQijdtap = (Uvozi*((2*B+Bsh)*Uvozi-B*tap*Uvozj*cos(di-dj-ph)+G*tap*Uvozj*sin(di-dj-ph)))/tap^3;
    dQjidtap = -((Uvozi*Uvozj*(B*cos(di-dj-ph)+G*sin(di-dj-ph)))/tap^2);
    
    %% dG/dU
    GpoU(indeks_P1 , k) = - dPijdtap; %dGfrom / dtap
    if tip_vozlisc(voz1) == 2 %PQ type of node
        GpoU(indeks_P1 + 1 , k) = - dQijdtap;
    end
    GpoU(indeks_P2 , k) = - dPjidtap;%dGto / dtap
    if tip_vozlisc(voz2) == 2 %PQ type of node
        GpoU(indeks_P2 + 1 , k) = - dQjidtap;
    end
    
    %% dH/dU
    HpoU(k,k) = 0;%odvod vozliscne napetosti po tap-u je 0
    if stevilo_pha_sh
        if find(phase_shifters(:,1) == veja)
            HpoU(find(phase_shifters(:,1) == veja) + stevilo_tap_ch , k) = dPijdtap;%dLinePower / dtap - in case a phase shifter and tap changer are in the same line
        end
    end
    
    %% dH/dX
    HpoX(k,indeks_reg_voz) = 1; %dRegNodeVolt / dVoltAmp = unity
end

for k=1:stevilo_pha_sh % partial derivatives for phases
    veja = phase_shifters(k,1);%tap branch
    voz1 = results.branch(veja,F_BUS);
    voz2 = results.branch(veja,T_BUS);
    indeks_P1 = ind_node_eq(voz1); %active power and volt. amp index for f node
    indeks_P2 = ind_node_eq(voz2); %active power and volt. amp index for t node
    
    Uvozi = results.bus(voz1,VM);
    Uvozj = results.bus(voz2,VM);
    di = results.bus(voz1,VA)/180*pi;
    dj = results.bus(voz2,VA)/180*pi;
    Z = results.branch(veja,BR_R) + 1i*results.branch(veja,BR_X);
    Bsh = results.branch(veja,BR_B);
    G = real(1/Z);
    B = imag(1/Z);
    tap = results.branch(veja,TAP);
    ph = results.branch(veja,SHIFT)/180*pi;

    dPijdph = (Uvozi*Uvozj*(B*cos(di-dj-ph)-G*sin(di-dj-ph)))/tap;
    dPjidph = -((Uvozi*Uvozj*(B*cos(di-dj-ph)+G*sin(di-dj-ph)))/tap);
    dPijUi = (2*G*Uvozi-tap*Uvozj*(G*cos(di-dj-ph)+B*sin(di-dj-ph)))/tap^2;
    dPijUj = -((Uvozi*(G*cos(di-dj-ph)+B*sin(di-dj-ph)))/tap);
    dPijdi = (Uvozi*Uvozj*(-B*cos(di-dj-ph)+G*sin(di-dj-ph)))/tap;
    dPijdj = (Uvozi*Uvozj*(B*cos(di-dj-ph)-G*sin(di-dj-ph)))/tap;
    dQijdph = (Uvozi*Uvozj*(G*cos(di-dj-ph)+B*sin(di-dj-ph)))/tap;
    dQjidph = (Uvozi*Uvozj*(-G*cos(di-dj-ph)+B*sin(di-dj-ph)))/tap;
    
    %% dG/dU
    GpoU(indeks_P1, k + stevilo_tap_ch) = - dPijdph;%dGfrom / dphase
    if tip_vozlisc(voz1) == 2 %PQ type of node
        GpoU(indeks_P1 + 1, k + stevilo_tap_ch) = - dQijdph;%dGfrom / dphase
    end
    GpoU(indeks_P2, k + stevilo_tap_ch) = - dPjidph;%dGto / dphase
    if tip_vozlisc(voz2) == 2 %PQ type of node
        GpoU(indeks_P2 + 1, k + stevilo_tap_ch) = - dQjidph;%dGto / dphase
    end
    
    %% dH/dU
    HpoU(k + stevilo_tap_ch, k + stevilo_tap_ch) = dPijdph;%dLinePower / dphase
    
    %% dH/dX
    if tip_vozlisc(voz1) == 2 %from bus is PQ type
        HpoX(k + stevilo_tap_ch, indeks_P1) = dPijUi;%dLinePower / dVoltAmp
        HpoX(k + stevilo_tap_ch, indeks_P1+1) = dPijdi;%dLinePower / dVoltAng
    else %from bus is PV type
        HpoX(k + stevilo_tap_ch, indeks_P1) = dPijdi;%dLinePower / dVoltAng
    end
    
    if tip_vozlisc(voz2) == 2 %same thing as abocve for to bus
        HpoX(k + stevilo_tap_ch, indeks_P2) = dPijUj;
        HpoX(k + stevilo_tap_ch, indeks_P2+1) = dPijdj;
    else
        HpoX(k + stevilo_tap_ch, indeks_P2) = dPijdj;
    end
    
end

sensitivity_of_H_to_U_taps_and_phasess = HpoU - HpoX * ( (GpoX) \ (GpoU) ); %total sensitivity of line active powers and nodal voltages w.r.t. all taps and phases
en