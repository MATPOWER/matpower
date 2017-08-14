function [g, dg] = opf_nle_fcn1(x)
Pg = x{1};
g = Pg(1)*Pg(2) - Pg(3);
if nargout == 2
    dg = sparse(1, length(Pg));
    dg(1:3) = [Pg(2) Pg(1) -1];
end
