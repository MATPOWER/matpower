function [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost
%IDX_COST   Defines constants for named column indices to gencost matrix.
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
%    4  NCOST       number of cost coefficients to follow for polynomial cost
%                   function, or number of data points for piecewise linear
%    5  COST        parameters defining total cost function begin in this col
%                  (MODEL = 1) : p0, f0, p1, f1, ..., pn, fn
%                       where p0 < p1 < ... < pn and the cost f(p) is defined
%                       by the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of
%                       the end/break-points of the piecewise linear cost
%                  (MODEL = 2) : cn, ..., c1, c0
%                       n+1 coefficients of an n-th order polynomial cost fcn,
%                       starting with highest order, where cost is
%                       f(p) = cn*p^2 + ... + c1*p + c0
% 
%   additional constants, used to assign/compare values in the MODEL column
%    1  PW_LINEAR   piecewise linear generator cost model
%    2  POLYNOMIAL  polynomial generator cost model
%
%   See also DEFINE_CONSTANTS.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
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

%% define cost models
PW_LINEAR   = 1;
POLYNOMIAL  = 2;

%% define the indices
MODEL       = 1;    %% cost model, 1 = piecewise linear, 2 = polynomial 
STARTUP     = 2;    %% startup cost in US dollars
SHUTDOWN    = 3;    %% shutdown cost in US dollars
NCOST       = 4;    %% number breakpoints in piecewise linear cost function,
                    %% or number of coefficients in polynomial cost function
COST        = 5;    %% parameters defining total cost function begin in this col
                    %% (MODEL = 1) : p0, f0, p1, f1, ..., pn, fn
                    %%      where p0 < p1 < ... < pn and the cost f(p) is defined
                    %%      by the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of
                    %%      the end/break-points of the piecewise linear cost
                    %% (MODEL = 2) : cn, ..., c1, c0
                    %%      n+1 coefficients of an n-th order polynomial cost fcn,
                    %%      starting with highest order, where cost is
                    %%      f(p) = cn*p^2 + ... + c1*p + c0
