classdef (Abstract) nme_wrapper_ac_nln < handle

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nme = [];           %% wrapped nm_element object
    end

    methods
        function obj = nme_wrapper_ac_nln_init(obj)
            obj.nme = obj.nme_class();      %% construct wrapped class
        end

        function build_nln_params(obj, nm, dm)
            %% build params for wrapped object
            obj.nme.build_params(nm, dm);

            %% remove other params
            obj.Y = [];
            obj.L = [];
            obj.M = [];
            obj.N = [];
            obj.i = [];
            obj.s = [];

            %% add nonlinear function/hessian params
            obj.inln = @(x_, sysx, idx)port_inj_current(obj.nme, x_, sysx, idx);
            obj.snln = @(x_, sysx, idx)port_inj_power(obj.nme, x_, sysx, idx);
            obj.inln_hess = @(x_, lam, sysx, idx)port_inj_current_hess(obj.nme, x_, lam, sysx, idx);
            obj.snln_hess = @(x_, lam, sysx, idx)port_inj_power_hess(obj.nme, x_, lam, sysx, idx);
        end
    end     %% methods
end         %% classdef
