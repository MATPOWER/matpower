function [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost
%IDX_COST   Defines variables for column indices to gencost.
%   [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% define cost models
PW_LINEAR	= 1;
POLYNOMIAL	= 2;

%% define the indices
MODEL		= 1;	%% cost model, 1 - piecewise linear, 2 - polynomial 
STARTUP		= 2;	%% startup cost in US dollars
SHUTDOWN	= 3;	%% shutdown cost in US dollars
N			= 4;	%% number breakpoints in piecewise linear cost function,
					%% or number of coefficients in polynomial cost function
COST		= 5;	%% beginning of cost parameters,
					%% piecewise linear data as:
					%%		x0, y0, x1, y1, x2, y2, ...
					%% and polynomial data as:
					%%		c2, c1, c0
					%% where the polynomial is c2*P^2 + c1*P + c0
					%% note: polynomials can be of any degree, highest
					%% order coefficient always goes first

return;
