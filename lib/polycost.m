function f = polycost(gencost, Pg, der)
%POLYCOST  Evaluates polynomial generator cost & derivatives
%   f = polycost(gencost, Pg)
%   df = polycost(gencost, Pg, 1)
%   d2f = polycost(gencost, Pg, 2)
%
%   gencost must contain only polynomial costs
%   Pg is in MW, not p.u. (works for Qg too)

%%----- initialize -----
%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

if nargin < 3
    der = 0;
end

if any(gencost(:, MODEL) == PW_LINEAR)
    error('polycost: all costs must be polynomial');
end

ng = length(Pg);
maxN = max(gencost(:, NCOST));
minN = min(gencost(:, NCOST));

%% form coefficient matrix where 1st column is constant term, 2nd linear, etc.
c = zeros(ng, maxN);
for n = minN:maxN
    k = find(gencost(:, NCOST) == n);   %% cost with n coefficients
    c(k, 1:n) = gencost(k, (COST+n-1):-1:COST);
end

%% do derivatives
for d = 1:der
    if size(c, 2) >= 2
        c = c(:, 2:maxN-d+1);
    else
        c = zeros(ng, 1);
        break;
    end
    for k = 2:maxN-d
        c(:, k) = k * c(:, k);
    end
end

%% evaluate polynomial
if isempty(c)
    f = zeros(size(Pg));
else
    f = c(:, 1);        %% constant term
    for k = 2:size(c, 2)
        f = f + c(:, k) .* Pg .^ (k-1);
    end
end
