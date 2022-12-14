classdef mme_line3p < mp.mm_element

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'line3p';
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            pp = nm.get_idx('port');
            nn = nme.np/2;

            %% branch active power flow
            for p = 1:nn
                s_fr = nm.soln.gs_(pp.i1.line3p(p):pp.iN.line3p(p)) * dm.base_kva;
                s_to = nm.soln.gs_(pp.i1.line3p(p+nn):pp.iN.line3p(p+nn)) * dm.base_kva;

                %% update in the data model
                dme.tab.(sprintf('pl%d_fr', p))(dme.on) = real(s_fr);
                dme.tab.(sprintf('ql%d_fr', p))(dme.on) = imag(s_fr);
                dme.tab.(sprintf('pl%d_to', p))(dme.on) = real(s_to);
                dme.tab.(sprintf('ql%d_to', p))(dme.on) = imag(s_to);
            end
        end
    end     %% methods
end         %% classdef
