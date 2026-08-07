function [results, stevilo_iteracij, success] = taps_and_phases_analysis(casedata, tap_changers_data, phase_shifters_data)
%taps_and_phases_analysis() caltulates taps and phases by iterating pwoer flows
%   [RESULTS, ITERATION_COUNT,SUCCESS] = taps_and_phases_analysis(casedata, tap_changers_data, phase_shifters_data)
%
%   Runs a Newtonian routine that calcualtes the taps and phases from 
%   specified nodal voltages and active power flow criteria. Taps and
%   phases not specified in the data specifically are fixed. Initial taps 
%   and phases are read from data. 
%
%   Inputs (all are mandatory), if there are no transformers use []:
%       casedata : A string containing the name of the file with the case
%       data (e.g. input = 'case300') 
%       tap_changers_data : Data of all the tap changers connected into the
%       system in the form of a matrix, the first column contains branches
%       the second contains regulated nodes and the third the voltage
%       magnitude in p.u.. For every transformer it is a line in the form 
%       of tap_changers_data = [branch , controlled node , votlage in p.u.]
%
%       phase_shifters_data : Data of all the phase shifters connected into
%       the system. The first column contains branches and the second
%       contains the active power flow into the line at the first node of
%       the branch of the device's insertion in p.u.. Every line is in the
%       form of phase_shifters_data = [branch , active power in p.u.]
%
%   Outputs (all are optional):
%       RESULTS : results struct 
%       ITERATION_COUNT : number of successive load flow calcualtions ran
%       which is also the iteration count with this method
%       SUCCESS : 0 if failed to converge, 1 if converged

%   Example:
%       results = taps_and_phases_analysis('case300',[1,9001,.95],[100,0]);
%   By Gorazd Bone, Faculty of Electrical Engineerinc, Ljubljana 


st_tap_ch = size(tap_changers_data,1); %number of tap changers
st_pha_sh = size(phase_shifters_data,1); %number of phase shifters

% U  %the control inputs are defined from initial conditions
H = zeros(st_tap_ch + st_pha_sh,1); %the criteria vector
itandp  = zeros(st_tap_ch + st_pha_sh,1); %initial taps and phases

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;   
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen; %#ok<*NASGU,*ASGLU>


%% Read data
mpc_ext = loadcase(casedata);
mpc_int = ext2int(mpc_ext);

%% Index reordering for tap changer control bus
ind_voz = [(1:size(mpc_ext.bus(:,BUS_I),1))' , mpc_ext.bus(:,BUS_I)];
for k=1: st_tap_ch
    tap_changers_data(k,2) = ind_voz(ind_voz(:,2) == tap_changers_data(k,2), 1);
end

%% Initial state determination
mpopt  =  mpoption('out.all',  0,  'verbose',  0);
results = ext2int(runpf(mpc_int,mpopt));
results.branch(results.branch(:,TAP)==0,TAP) = 1; % 0 initial tap is regarded as 1
Sbazni = mpc_int.baseMVA; %razmerje kako so podane moci v vektorju results.branch
if size(tap_changers_data,1>0)
    itandp(1 : st_tap_ch) = results.branch(tap_changers_data(:,1),TAP);
    H(1 : st_tap_ch) = results.bus(tap_changers_data(:,2),VM) - tap_changers_data(:,3);
end
if size(phase_shifters_data,1>0)
    itandp(st_tap_ch + 1 : st_tap_ch + st_pha_sh) = results.branch(phase_shifters_data(:,1),SHIFT)/180*pi;
    H(st_tap_ch + 1 : st_tap_ch + st_pha_sh) = results.branch(phase_shifters_data(:,1),PF)/Sbazni - phase_shifters_data(:,2);
end
U = itandp; % the input vector

%% Read node types
pqtip = results.bus(:,BUS_TYPE)==PQ;
pvtip = results.bus(:,BUS_TYPE)==PV;
tip_vozlisc = pqtip*2 + pvtip; % nodal type, 2 for PQ 1 for PV, corresponds to the nubmer of equations

%% Iterate taps and phases
stevilo_iteracij = 0; %iteration count
while (norm(H,2)>1e-5) && (stevilo_iteracij<40)
    J = shiftJac_taps_phases(results); % jacobian matrix
    obcutljivosti = sensit_to_taps_and_phases( results, tap_changers_data, phase_shifters_data, J , tip_vozlisc);% system sensitivity calculation
    U = U - obcutljivosti\H;%step by Newton's method
    
    if st_tap_ch
        results.branch(tap_changers_data(:,1),TAP) = U(1:st_tap_ch); %apply new taps
    end
    if st_pha_sh
        results.branch(phase_shifters_data(:,1),SHIFT) = U(st_tap_ch + 1 : st_tap_ch + st_pha_sh)*180/pi;%apply new phases
    end
    
    results = ext2int(runpf(results,mpopt)); % rerun load flow with new values
    if size(tap_changers_data,1>0)
        H(1 : st_tap_ch) = results.bus(tap_changers_data(:,2),VM) - tap_changers_data(:,3); % criteria for tap changers
    end
    if size(phase_shifters_data,1>0)
        H(st_tap_ch + 1 : st_tap_ch + st_pha_sh) = results.branch(phase_shifters_data(:,1),PF)/Sbazni - phase_shifters_data(:,2);% criteria for phase shifters
    end
    stevilo_iteracij = stevilo_iteracij + 1;
end

results = int2ext(results);

success = (~(norm(H,2)>1e-5));