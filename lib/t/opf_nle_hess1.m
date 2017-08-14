function d2G = opf_nle_hess1(x, lambda)
Pg = x{1};
n = length(Pg);
d2G = sparse([1;2], [2;1], [lambda; lambda], n, n);
end
