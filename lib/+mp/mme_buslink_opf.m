classdef (Abstract) mme_buslink_opf < mp.mme_buslink
% mp.mme_buslink_opf - Math model element abstract base class for 1-to-3-phase buslink for OPF.
%
% Abstract math model element base class for 1-to-3-phase buslink elements
% for OPF problems.
%
% Implements (currently empty) method for forming an interior initial point.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %
        end
        
        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            ss = nm.get_idx('state');

            for p = 1:nme.nz
                %% buslink complex power flows
                zbl = nm.soln.z(ss.i1.buslink(p):ss.iN.buslink(p));                

                %% update in the data model
                dme.tab.(sprintf('pg%d_start', p))(dme.on) = real(zbl) * dm.base_kva;
                dme.tab.(sprintf('qg%d_start', p))(dme.on) = imag(zbl) * dm.base_kva;
            end
        end
    end     %% methods
end         %% classdef
