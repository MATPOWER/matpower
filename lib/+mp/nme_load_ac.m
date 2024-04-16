classdef (Abstract) nme_load_ac < mp.nme_load% & mp.form_ac
% mp.nme_load_ac - Network model element abstract base class for load for AC formulations.
%
% Builds the parameters :math:`\sv` and :math:`\YY` and nonlinear functions
% :math:`\Snln(\X)` and :math:`\Inln(\X)`.

%   MATPOWER
%   Copyright (c) 2019-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build_params(obj, nm, dm)
            %
            build_params@mp.nme_load(obj, nm, dm);  %% call parent

            dme = obj.data_model_element(dm);

            %% constant complex power demand
            obj.s = dme.pd + 1j * dme.qd;

            %% nominal complex power from constant current demand
            %% (scaled by voltage magnitude)
            if any(dme.pd_i) || any(dme.qd_i)
%                 obj.i = pd_i - 1j qd_i;   %% power as function of complex voltage, not voltage magnitude (as desired)
                Sd = dme.pd_i + 1j * dme.qd_i;
                obj.snln = @(x_, sysx, idx)port_inj_power_nln(obj, Sd, x_, sysx, idx);
                obj.inln = @(x_, sysx, idx)port_inj_current_nln(obj, Sd, x_, sysx, idx);
            end

            %% nominal complex power from constant impedance demand
            %% (scaled by voltage magnitude squared)
            if any(dme.pd_z) || any(dme.qd_z)
                Y = dme.pd_z - 1j * dme.qd_z;
                nd = length(Y);
                obj.Y = sparse(1:nd, 1:nd, Y, nd, nd);
            end
        end

        function [I, Iv1, Iv2, Izr, Izi] = port_inj_current_nln(obj, Sd, x_, sysx, idx)
            %
            error('mp.nme_load_ac.port_inj_current_nln: not yet implemented');
%             if nargin < 5
%                 idx = [];
%                 if nargin < 4
%                     sysx = 1;
%                 end
%             end
% 
%             [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);
%             if isempty(idx)
%                 Sdi = Sd;
%             else
%                 Sdi = Sd(idx);
%             end
%             S = abs(vi_) .* Sdi;
%             I = S ./ vi_;
% 
%             if nargout > 1
%                 nv = length(v_);
%                 nz = length(z_);
%                 ni = length(S);
%                 if isempty(idx)
%                     idx = (1:ni);
%                 end
%                 Sv1 = sparse(ni, nv);
%                 Sv2 = sparse(1:ni, idx, Sdi, ni, nv);
%                 if nargout > 3
%                     Szr = sparse(ni, nz);
%                     Szi = Szr;
%                 end
%                 if sysx
%                     Ct = obj.C';
%                     Sv1 = Sv1 * Ct;
%                     Sv2 = Sv2 * Ct;
%                     if nargout > 3  %% Szr, Szi are empty, but num of rows is needed
%                         Dt = obj.D';
%                         Szr = Szr * Dt;
%                         Szi = Szi * Dt;
%                     end
%                 end
%             end
        end

        function [S, Sv1, Sv2, Szr, Szi] = port_inj_power_nln(obj, Sd, x_, sysx, idx)
            %
            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);
            if isempty(idx)
                Sdi = Sd;
            else
                Sdi = Sd(idx);
            end
            S = abs(vi_) .* Sdi;

            if nargout > 1
                nv = length(v_);
                nz = length(z_);
                ni = length(S);
                if isempty(idx)
                    idx = (1:ni);
                end
                Sv1 = sparse(ni, nv);
                Sv2 = sparse(1:ni, idx, Sdi, ni, nv);
                if nargout > 3
                    Szr = sparse(ni, nz);
                    Szi = Szr;
                end
                if sysx
                    Ct = obj.C';
                    Sv1 = Sv1 * Ct;
                    Sv2 = Sv2 * Ct;
                    if nargout > 3  %% Szr, Szi are empty, but num of rows is needed
                        Dt = obj.D';
                        Szr = Szr * Dt;
                        Szi = Szi * Dt;
                    end
                end
            end
        end
    end     %% methods
end         %% classdef
