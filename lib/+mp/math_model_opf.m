classdef (Abstract) math_model_opf < mp.math_model
% mp.math_model_opf - Abstract base class for optimal power flow (OPF) **math model** objects.
%
% Provide implementations for adding system variables to the mathematical
% model and creating an interior starting point.

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
        function tag = task_tag(obj)
            %

            tag = 'opf';
        end

        function name = task_name(obj)
            %

            name = 'Optimal Power Flow';
        end

        function ad = build_aux_data(obj, nm, dm, mpopt)
            %

            %% create aux_data struct
            ad = obj.build_base_aux_data(nm, dm, mpopt);
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
            %

            %% add network voltage and non-voltage state variables
            vars = horzcat(nm.model_vvars(), nm.model_zvars());
            for vtype = vars
                st = nm.(vtype{1});     %% set type
                d = st.data;
                mmx_i1 = obj.var.N + 1;
                for k = 1:st.NS
                    name = st.order(k).name;
                    idx = st.order(k).idx;
                    if isempty(idx)
                        obj.add_var(name, st.idx.N.(name), d.v0.(name), d.vl.(name), d.vu.(name), d.vt.(name));
                    else
                        if all(cell2mat(idx) == 1)
                            dim = size(st.idx.N.(name));
                            if dim(end) == 1, dim(end) = []; end   %% delete trailing 1
                            obj.init_indexed_name('var', name, num2cell(dim));
                        end
                        sn = struct('type', {'()'}, 'subs', {idx});
                        sc = sn;
                        sc.type = '{}';
                        N = subsref(st.idx.N.(name), sn);
                        v0 = subsref(d.v0.(name), sc);
                        vl = subsref(d.vl.(name), sc);
                        vu = subsref(d.vu.(name), sc);
                        vt = subsref(d.vt.(name), sc);
                        obj.add_var(name, idx, N, v0, vl, vu, vt);
                    end
                end
                mmx_iN = obj.var.N;
                obj.aux_data.var_map{end+1} = ...
                    {vtype{1}, [], [], [], mmx_i1, mmx_iN, []};
            end
        end

        function x0 = interior_x0(obj, mm, nm, dm, x0)
            %

            if nargin < 5 || isempty(x0)
                %% generic interior point
                [x0, xmin, xmax] = mm.params_var();     %% init var & bounds
                s = 1;                      %% set init point inside bounds by s
                lb = xmin; ub = xmax;
                lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
                ub(xmax ==  Inf) =  1e10;   %% temporarily to avoid errors in next line
                x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
                k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
                x0(k) = xmax(k) - s;                    %% set just below upper bound
                k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
                x0(k) = xmin(k) + s;                    %% set just above lower bound
            end

            %% allow each node or state creating element to initialize its vars
            for k = 1:length(nm.elements)
                nme = nm.elements{k};   %% network model element
                if nme.nn || nme.nz     %% element creates nodes or states
                    mme = nme.math_model_element(mm);   %% math model element
                    if ~isempty(mme)
                        x0 = mme.interior_x0(mm, nm, dm, x0);
                    end
                end
            end
        end

        function varef1 = interior_va(obj, nm, dm)
            %

            %% return scalar va equal to angle of first reference node
            ad = obj.aux_data;
            ref1 = ad.ref(1);
            varef1 = ad.va(ref1);
        end
    end     %% methods
end         %% classdef
