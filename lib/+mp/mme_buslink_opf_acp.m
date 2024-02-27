classdef mme_buslink_opf_acp < mp.mme_buslink_opf
% mp.mme_buslink_opf_acp - Math model element for 1-to-3-phase buslink for AC polar voltage OPF.
%
% Math model element class for 1-to-3-phase buslink elements for AC polar
% OPF problems.
%
% Implements method for adding constraints to match voltages across each
% buslink.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %
            nme = obj.network_model_element(nm);

            %% voltage equality constraints
            [A, b_va, b_vm] = nme.voltage_constraints();
            idx = {{}, {1}, {2}, {3}};
            vs_va = struct('name', {'Va', 'Va3', 'Va3', 'Va3'}, 'idx', idx);
            vs_vm = struct('name', {'Vm', 'Vm3', 'Vm3', 'Vm3'}, 'idx', idx);
            mm.add_lin_constraint('buslink_va', A, b_va, b_va, vs_va);
            mm.add_lin_constraint('buslink_vm', A, b_vm, b_vm, vs_vm);
        end
    end     %% methods
end         %% classdef
