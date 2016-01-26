function [mpco, NewGenIdx] = addgen2mpc(mpci, gen, gencost, fuel)
%ADDGEN2MPC Adds a set of generators to a MATPOWER case struct.
%
%   [NEW_MPC, IDX] = ADDGEN2MPC(MPC, GEN, GENCOST, FUEL)
%
%   Inserts a set of generators to a MATPOWER case struct by adding the
%   info of the new generators (in inputs GEN, GENCOST and FUEL) to the
%   bottom of the corresponding tables (fields of MPC). Dimensions must
%   be consistent.
%
%   Inputs:
%       MPC : standard MATPOWER case struct, with the following additional
%             fields:
%               .genfuel : (optional) cell array of strings with fuel type
%                          (filled with 'unknown' if missing)
%               .i<fuel> : vector of generator indices for units of the
%                          specified fuel type
%       GEN: standard MATPOWER generator matrix for units to be added
%       GENCOST: standard MATPOWER generator cost matrix for units to be added
%       FUEL: string or cell array of names of fuel or type of new generators
%
%   Outputs:
%       NEW_MPC : the new MATPOWER case with the additional generators
%                 appended to GEN, GENCOST, GENFUEL and additional field
%                   .i<fuel> : vector of generator indices for units of the
%                              specified fuel type
%       IDX : generator indices of the newly added generators

%   MPSOPF
%   $Id$
%   by Daniel Munoz-Alvarez and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2013 by Power System Engineering Research Center (PSERC)

[ng, ncg]     = size(mpci.gen);
[nr, nc]      = size(gen);
[ng2, ncgc]   = size(mpci.gencost);
[nr2, ncgc2]  = size(gencost);
if nr ~= nr2
    error('addgen2mpc: number of rows in GEN (%d) and GENCOST (%d) must agree', nr, nr2);
end
if ng ~= ng2
    error('addgen2mpc: number of rows in MPC.GEN (%d) and MPC.GENCOST (%d) must agree', ng, ng2);
end

%% append rows to GEN
mpco = mpci;
if nc <= ncg
    mpco.gen = [mpco.gen ; gen zeros(nr, ncg-nc)];
else
    error('addgen2mpc: number of columns of GEN (%d) exceeds number of columns of MPC.GEN (%d)', nc, ncg);
end

%% append rows to GENCOST
dim_gencost = size(gencost);
dim_gencost_cur = size(mpco.gencost);
if dim_gencost(2) < dim_gencost_cur(2)
    mpco.gencost = [mpco.gencost ; [gencost zeros(dim_gencost(1),dim_gencost_cur(2)-dim_gencost(2))]];
elseif dim_gencost(2) > dim_gencost_cur(2)
    mpco.gencost = [mpco.gencost zeros(dim_gencost_cur(1),dim_gencost(2)-dim_gencost_cur(2))];
    mpco.gencost = [mpco.gencost ; gencost];
else
    mpco.gencost = [mpco.gencost ; gencost];
end

%% initialize 'genfuel'  and 'i<fuel>' fields if missing
if ~isfield(mpco, 'genfuel')
    mpco.genfuel = mat2cell(repmat('unknown', ng, 1), ones(ng,1), 7);
end
if ~isfield(mpco,['i' fuel])
    mpco.(['i' fuel]) = [];
end

%% append 'genfuel' and 'i<fuel>' fields
if iscell(fuel)
    for i = dim_gencost(1)
        mpco.genfuel = [mpco.genfuel ; fuel{i}];
        mpco.(['i' fuel{i}]) = [ mpco.(['i' fuel{i}]) ; dim_gencost_cur(1) + i ];
    end
elseif ischar(fuel)
    for i = 1:dim_gencost(1)
        mpco.genfuel = [ mpco.genfuel ; fuel ];
        mpco.(['i' fuel]) = [ mpco.(['i' fuel]) ; dim_gencost_cur(1) + i ];
    end
else
    error('addgen2mpc: FUEL must be a string (or cell array of strings) indicating the fuel type of the new generators');
end

NewGenIdx = ( ng + 1 : size(mpco.gen,1) )';
