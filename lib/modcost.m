function gencost = modcost(gencost, alpha, modtype)
%MODCOST  Modifies generator costs by shifting or scaling (F or X).
%   NEWGENCOST = MODCOST(GENCOST, ALPHA)
%   NEWGENCOST = MODCOST(GENCOST, ALPHA, MODTYPE)
%
%   For each generator cost F(X) (for real or reactive power) in
%   GENCOST, this function modifies the cost by scaling or shifting
%   the function by ALPHA, depending on the value of MODTYPE, and
%   and returns the modified GENCOST. Rows of GENCOST can be a mix
%   of polynomial or piecewise linear costs. ALPHA can be a scalar,
%   applied to each row of GENCOST, or an NG x 1 vector, where each
%   element is applied to the corresponding row of GENCOST.
%
%   MODTYPE takes one of the 4 possible values (let F_alpha(X) denote the
%   the modified function):
%       'SCALE_F' (default) : F_alpha(X)         == F(X) * ALPHA
%       'SCALE_X'           : F_alpha(X * ALPHA) == F(X)
%       'SHIFT_F'           : F_alpha(X)         == F(X) + ALPHA
%       'SHIFT_X'           : F_alpha(X + ALPHA) == F(X)

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define named indices into data matrices
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

if nargin < 3
    modtype = 'SCALE_F';
end

[ng, m] = size(gencost);

if ng ~= 0
    if length(alpha) ~= ng
        if length(alpha) == 1 && ng > 1     %% scalar, make it a col vector
            alpha = alpha * ones(ng, 1);
        else
            error('modcost: ALPHA must be a scalar or col vector with NG rows');
        end
    elseif size(alpha, 2) ~= 1
        alpha = alpha';                     %% convert row vector to col vector
    end

    ipwl = find(gencost(:, MODEL) == PW_LINEAR);
    ipol = find(gencost(:, MODEL) == POLYNOMIAL);
    c = gencost(ipol, COST:m);

    switch modtype
        case 'SCALE_F',
            gencost(ipol, COST:m)       = diag(alpha(ipol)) * c;
            gencost(ipwl, COST+1:2:m)   = diag(alpha(ipwl)) * gencost(ipwl, COST+1:2:m);
        case 'SCALE_X',
            for k = 1:length(ipol)
                n = gencost(ipol(k), NCOST);
                for i = 1:n
                    gencost(ipol(k), COST+i-1) = c(k, i) / alpha(ipol(k))^(n-i);
                end
            end
            gencost(ipwl, COST:2:m-1)   = diag(alpha(ipwl)) * gencost(ipwl, COST:2:m-1);
        case 'SHIFT_F',
            for k = 1:length(ipol)
                n = gencost(ipol(k), NCOST);
                gencost(ipol(k), COST+n-1) = alpha(ipol(k)) + c(k, n);
            end
            gencost(ipwl, COST+1:2:m)   = ...
                diag(alpha(ipwl)) * ones(length(ipwl), (m+1-COST)/2) + ...
                    gencost(ipwl, COST+1:2:m);
        case 'SHIFT_X',
            for k = 1:length(ipol)
                n = gencost(ipol(k), NCOST);
                gencost(ipol(k), COST:COST+n-1) = polyshift(c(k, 1:n)', alpha(ipol(k)))';
            end
            gencost(ipwl, COST:2:m-1)   = ...
                diag(alpha(ipwl)) * ones(length(ipwl), (m+1-COST)/2) + ...
                    gencost(ipwl, COST:2:m-1);
        otherwise
            error('modcost: ''%s'' is not a valid modtype\n', modtype);
    end
end


%%-----  POLYSHIFT  -----
function d = polyshift(c, a)
%POLYSHIFT  Returns the coefficients of a horizontally shifted polynomial.
%
%   D = POLYSHIFT(C, A) shifts to the right by A, the polynomial whose
%   coefficients are given in the column vector C.
%
%   Example: For any polynomial with n coefficients in c, and any values
%   for x and shift a, the f - f0 should be zero.
%       x = rand;
%       a = rand;
%       c = rand(n, 1);
%       f0 = polyval(c, x)
%       f  = polyval(polyshift(c, a), x+a)

n = length(c);
d = zeros(size(c));
A = (-a * ones(n, 1)) .^ ((0:n-1)');
b = ones(n, 1);
for k = 1:n
    d(n-k+1) = b' * ( c(n-k+1:-1:1) .* A(1:n-k+1) );
    b = cumsum(b(1:n-k));
end
