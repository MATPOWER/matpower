function [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost
% idx_cost - Defines constants for named column indices to ``gencost`` matrix.
% ::
%
%   Example:
%
%   [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    start = gencost(4, STARTUP);       % get startup cost of generator 4
%    gencost(2, [MODEL, NCOST:COST+1]) = [ POLYNOMIAL 2 30 0 ];
%    % set the cost of generator 2 to a linear function COST = 30 * Pg
% 
%   The index, name and meaning of each column of the gencost matrix is given
%   below:
% 
%   columns 1-5
%    1  MODEL       cost model, 1 = piecewise linear, 2 = polynomial
%    2  STARTUP     startup cost in US dollars
%    3  SHUTDOWN    shutdown cost in US dollars
%    4  NCOST       number N = n+1 of data points to follow defining an
%                   n-segment piecewise linear cost function, or of cost
%                   coefficients defining an n-th order polynomial cost function
%    5  COST        parameters defining total cost function f(p) begin in
%                   this column
%                   (MODEL = 1) : p1, f1, p2, f2, ..., pN, fN
%                       where p1 < p2 < ... < pN and the cost f(p) is defined
%                       by the coordinates (p1,f1), (p2,f2), ..., (pN,fN) of
%                       the end/break-points of the piecewise linear cost fcn
%                   (MODEL = 2) : cn, ..., c1, c0
%                       N coefficients of an n-th order polynomial cost
%                       function, starting with highest order, where cost is
%                       f(p) = cn*p^n + ... + c1*p + c0
%
%   additional constants, used to assign/compare values in the MODEL column
%    1  PW_LINEAR   piecewise linear generator cost model
%    2  POLYNOMIAL  polynomial generator cost model
%
% See also define_constants.

%   MATPOWER
%   Copyright (c) 1996-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% define cost models
PW_LINEAR   = 1;
POLYNOMIAL  = 2;

%% define the indices
MODEL       = 1;    %% cost model, 1 = piecewise linear, 2 = polynomial 
STARTUP     = 2;    %% startup cost in US dollars
SHUTDOWN    = 3;    %% shutdown cost in US dollars
NCOST       = 4;    %% number N = n+1 of end/breakpoints in piecewise linear
                    %% cost function, or of coefficients in polynomial cost fcn
COST        = 5;    %% parameters defining total cost function begin in this col
                    %% (MODEL = 1) : p1, f1, p2, f2, ..., pN, fN
                    %%      where p1 < p2 < ... < pN and the cost f(p) is defined
                    %%      by the coordinates (p1,f1), (p2,f2), ..., (pN,fN) of
                    %%      the end/break-points of the piecewise linear cost
                    %% (MODEL = 2) : cn, ..., c1, c0
                    %%      N coefficients of an n-th order polynomial cost fcn,
                    %%      starting with highest order, where cost is
                    %%      f(p) = cn*p^n + ... + c1*p + c0
