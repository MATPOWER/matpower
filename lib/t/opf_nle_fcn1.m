function [g, dg] = opf_nle_fcn1(x)
% opf_nle_fcn1 - Example user-defined nonlinear OPF constraint function.
Pg = x{1};
g = Pg(1)*Pg(2) - Pg(6);
if nargout == 2
    dg = [Pg(2) Pg(1) 0 0 0 -1];
end
