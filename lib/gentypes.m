function gt = gentypes()
%GENTYPES   Return list of standard values for generator unit types
%
%   GT = GENTYPES()
%
%   Returns a cell array of strings containing the following standard
%   two character generator unit types for use in the optional MPC.GENTYPE
%   field of the MATPOWER case struct. This is to be considered an unordered
%   list, where the position of a particular fuel type in the list is not
%   defined and is therefore subject to change.
%
%   From Form EIA-860 Instructions, Table 2. Prime Mover Codes and Descriptions
%   https://www.eia.gov/survey/form/eia_860/instructions.pdf
%       BA  - Energy Storage, Battery
%       CE  - Energy Storage, Compressed Air
%       CP  - Energy Storage, Concentrated Solar Power
%       FW  - Energy Storage, Flywheel
%       PS  - Hydraulic Turbine, Reversible (pumped storage)
%       ES  - Energy Storage, Other
%       ST  - Steam Turbine, including nuclear, geothermal and solar steam
%               (does not include combined cycle)
%       GT  - Combustion (Gas) Turbine (includes jet engine design)
%       IC  - Internal Combustion Engine (diesel, piston, reciprocating)
%       CA  - Combined Cycle Steam Part
%       CT  - Combined Cycle Combustion Turbine Part
%               (type of coal or solid must be reported as energy source
%               for integrated coal gasification)
%       CS  - Combined Cycle Single Shaft
%               (combustion turbine and steam turbine share a single generator)
%       CC  - Combined Cycle Total Unit
%               (use only for plants/generators that are in planning stage,
%               for which specific generator details cannot be provided)
%       HA  - Hydrokinetic, Axial Flow Turbine
%       HB  - Hydrokinetic, Wave Buoy
%       HK  - Hydrokinetic, Other
%       HY  - Hydroelectric Turbine (includes turbines associated with
%               delivery of water by pipeline)
%       BT  - Turbines Used in a Binary Cycle
%               (including those used for geothermal applications)
%       PV  - Photovoltaic
%       WT  - Wind Turbine, Onshore
%       WS  - Wind Turbine, Offshore
%       FC  - Fuel Cell
%       OT  - Other
%   Additional codes (some from PowerWorld)
%       UN  - Unknown
%       JE  - Jet Engine
%       NB  - ST - Boiling Water Nuclear Reactor
%       NG  - ST - Graphite Nuclear Reactor
%       NH  - ST - High Temperature Gas Nuclear Reactor
%       NP  - ST - Pressurized Water Nuclear Reactor
%       IT  - Internal Combustion Turbo Charged
%       SC  - Synchronous Condenser
%       DC  - represents DC ties
%       MP  - Motor/Pump
%       W1  - Wind Turbine, Type 1
%       W2  - Wind Turbine, Type 2
%       W3  - Wind Turbine, Type 3
%       W4  - Wind Turbine, Type 4
%       SV  - Static Var Compensator
%       DL  - Dispatchable Load
%
%   Example:
%       if ~ismember(mpc.gentype{k}, gentypes())
%           error('unknown generator unit type');
%       end
%
%   See also GENFUELS, SAVECASE.

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

gt = { ...
    'BA';
    'CE';
    'CP';
    'FW';
    'PS';
    'ES';
    'ST';
    'GT';
    'IC';
    'CA';
    'CT';
    'CS';
    'CC';
    'HA';
    'HB';
    'HK';
    'HY';
    'BT';
    'PV';
    'WT';
    'WS';
    'FC';
    'OT';
    'UN';
    'JE';
    'NB';
    'NG';
    'NH';
    'NP';
    'IT';
    'SC';
    'DC';
    'MP';
    'W1';
    'W2';
    'W3';
    'W4';
    'SV';
    'DL';
};
