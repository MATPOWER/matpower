function Zpv = make_zpv(pv,nb,nl,f,Zb)
%MAKE_ZPV  Calculates loop impedances for all PV buses.
%
%   Zpv = make_zpv(pv,nb,nl,f,Zb)
%
%   Loop impedance of a PV bus is defined as impedance of the path between
%   the bus and the slack bus. The mutual impedance between two PV buses is
%   the impedance of the joint part of the two path going from each of the
%   PV buses to the slack bus. The impedances are calculated as bus voltages
%   in cases when at one of the PV buses we inject current of -1 A. All
%   voltages are calculated with the backward-forward sweep method. The
%   input variables are the vector of indicies with "from" buses for each
%   branch, the vector of branch impedances and indicies of PV buses.
%
%   See also CALC_V_PQ_SUM.

npv = length(pv);
Zpv = zeros(npv);
for ipv = 1:npv
    V = zeros(nb,1);
    I = zeros(nb,1);
    I(pv(ipv)) = -1;
    for k = nl:-1:2
        i = f(k);
        I(i) = I(i) + I(k);
    end
    for k = 2:nl
        i = f(k);
        V(k) = V(i) - Zb(k) * I(k);
    end 
    Zpv(:,ipv) = V(pv);
end