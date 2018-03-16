function d2VaDif = opf_branch_ang_hess(x, lambda, Aang, lang, uang, iang, mpopt)
%% unpack data
[Vr, Vi] = deal(x{:});

%% problem dimensions
nb = length(Vi);            %% number of buses

%% ----- evaluate constraint gradients -----
nlam = length(lambda) / 2;
if nlam
    lamUp = lambda(1:nlam);
    lamLow = lambda((1:nlam)+nlam);
else
    lamUp = zeros(0,1);   
    lamLow = zeros(0,1);  
end
%% ----- evaluate constraint gradients -----
diagVrVm4   = sparse(1:nb, 1:nb, 2*Vr./(Vr.^4 + 2*(Vr.^2).*(Vi.^2)+ Vi.^4), nb, nb);
diagViVm4   = sparse(1:nb, 1:nb, 2*Vi./(Vr.^4 + 2*(Vr.^2).*(Vi.^2)+ Vi.^4), nb, nb);
diagVm2     = sparse(1:nb, 1:nb, 1./(Vr.^2 + Vi.^2), nb, nb);
diagVr = sparse(1:nb, 1:nb, Vr, nb, nb);
diagVi = sparse(1:nb, 1:nb,Vi, nb, nb);
%% for upper limit
    diagAlam = sparse(1:nb, 1:nb, Aang'*lamUp, nb, nb);
    diagAlamV2 = sparse(1:nb, 1:nb, diagVm2*Aang'*lamUp, nb, nb);    
    
    VaDifU_ii = - diagVr*diagAlam*diagViVm4;
    VaDifU_ir = diagAlamV2 - diagVi*diagAlam*diagVrVm4;
    VaDifU_ri = - VaDifU_ir;
    VaDifU_rr = diagVi*diagAlam*diagVrVm4;
%% for lower limit
    diagAlam = sparse(1:nb, 1:nb, Aang'*lamLow, nb, nb);
    diagAlamV2 = sparse(1:nb, 1:nb, diagVm2*Aang'*lamLow, nb, nb);    
    
    VaDifL_ii = - diagVr*diagAlam*diagViVm4;
    VaDifL_ir = diagAlamV2 - diagVi*diagAlam*diagVrVm4;
    VaDifL_ri = - VaDifL_ir;
    VaDifL_rr = diagVi*diagAlam*diagVrVm4; 
%% construct Hessian
d2VaDif = - [ VaDifL_ii VaDifL_ir; VaDifL_ri VaDifL_rr] + [ VaDifU_ii VaDifU_ir; VaDifU_ri VaDifU_rr];        
end