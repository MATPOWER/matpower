function [g, dg] = opf_nle_fcn1(x)
Pg = x{1};
g = Pg(1)*Pg(2) - Pg(6);
if nargout == 2
    dg = [Pg(2) Pg(1) 0 0 0 -1];
end
