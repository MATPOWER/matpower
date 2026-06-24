function [ J ] = shiftJac_taps_phases( results )
%shiftJac_taps_phases  Constructs the Jacobian matrix using Matpower 
%supplied makeJac and reorderes it for taps and phases analysis.
%   [J ] = shiftJac_taps_phases( results )
%   Inputs :
%       results: results from a runpf(), interanly idnexed (use ext2int)
%
%   Output :
%       J : Jacobian matrix for power mismatch load flow, ordereing
%       appropriate fo sensit_to_taps_and_phases()
%
%   by Gorazd Bone, Faculty of Electrical Engineering, Ljubljana 

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;    
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...                    
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch; %#ok<*ASGLU,*NASGU>
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;  

pq = results.bus(:,BUS_TYPE)==PQ;
pv = results.bus(:,BUS_TYPE)==PV;
ref = find(pq + pv == 0);

tip_vozlisc = pq * 2 + pv;
ind_node_eq = tip_vozlisc * 0;
for k=1:size(tip_vozlisc,1)
    ind_node_eq(k) = sum(tip_vozlisc(1:k-1)) + 1;
end
loc_P_pv = ind_node_eq(pv==1);
loc_P_pq = ind_node_eq(pq==1);
loc_Q_pq = loc_P_pq+1;


per_rows = sparse([loc_P_pv;loc_P_pq;loc_Q_pq] , 1:sum(pq * 2 + pv),1); 
per_cols = sparse(1:sum(pq * 2 + pv) , [loc_P_pv;loc_Q_pq;loc_P_pq] , 1);
%
%   MatPower Jacobi has the following oredering for the missmatch:
%
%       P_vect_of_pvnodes    
%       P_vect_of_pqnodes    
%       Q_vect_of_pqnodes    
%
% and for the unknowns:
%
%       angle_vect_of_pvnodes    
%       angle_vect_of_pqnodes    
%       voltage_vect_of_pqnodes    
%
% the used Jacobi has the ordering node-by-node, for the missmatch
% 
%       P_of_node1
%       Q_of_node1 (if PQ)
%       P_of_node2
%       Q_of_node2 (if PQ)
%
% and for the unknowns:
%
%       voltage_of_node1 (if PQ)
%       anle_of_node1
%       voltage_of_node2 (if PQ)
%       anle_of_node2 (if PQ)
%
J = - per_rows * makeJac(results.baseMVA, results.bus, results.branch, results.gen) * per_cols;

