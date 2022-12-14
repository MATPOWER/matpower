classdef data_model < mp.element_container
%MP.DATA_MODEL  Base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        base_mva        %% system per unit MVA base [1]
        base_kva        %% system per unit kVA base [2]
        source          %% source of data (e.g. mpc, MATPOWER case struct)
        userdata = struct();
    end     %% properties
    %% [1]  for balanced single-phase systems/sections, must be provided if
    %%      system includes any 'bus' elements
    %% [2]  for unbalanced 3-phase systems/sections, must be provided if
    %%      system includes any 'bus3p' elements

    methods
        %% constructor
        function obj = data_model()
            %% call parent constructor
            obj@mp.element_container();
            obj.element_classes = ...
                { @mp.dme_bus, @mp.dme_gen, @mp.dme_load, ...
                    @mp.dme_branch, @mp.dme_shunt };
        end

        function new_obj = copy(obj)
            %% make shallow copy of object
            new_obj = eval(class(obj));  %% create new object
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end

            %% make copies of each individual element
            new_obj.elements = new_obj.elements.copy();
        end

        function obj = build(obj, d, dmc)
            %% create empty element objects for each class
            obj.elements = mp.mapped_array();
            for k = 1:length(obj.element_classes)
                dme_class = obj.element_classes{k};
                dme = dme_class();      %% element constructor
                obj.elements.add_elements(dme, dme.name);
            end

            %% import data from external format
            dmc.import(obj, d);

            %% count element objects for each class
            %% remove if count is zero
            for k = length(obj.elements):-1:1
                obj.elements{k}.count(obj);     %% set nr
                if obj.elements{k}.nr == 0
                    obj.elements.delete_elements(k);
                end
            end

            obj.initialize();       %% get IDs, self status
            obj.update_status();    %% update status wrt connectivity, define on/off
            obj.build_params();     %% get parameters
        end

        function obj = initialize(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.initialize(obj);
            end
        end

        function obj = update_status(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.update_status(obj);
            end
        end

        function obj = build_params(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.build_params(obj);
            end
        end

        function n = online(obj, name)
            if obj.elements.is_index_name(name)
                n = obj.elements.(name).n;
            else
                n = 0;
            end
        end

        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA MODEL CLASS : %s\n', class(obj));
            fprintf('        base_mva : %g\n', obj.base_mva);

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name            nr         n      class\n');
            fprintf('-- -----------   --------  --------  --------------------\n');
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                fprintf('%2d  %-13s %6d %9d    %s\n', k, dme.name, dme.nr, dme.n, class(dme));
%                 dme
            end

            %% user data
            fields = fieldnames(obj.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                         size       class\n');
                fprintf(' ------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = obj.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-24s  %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            end
        end

        function [out, add] = pp_flags(obj, mpopt)
            %% output control flags
            ne = length(obj.elements);
            names = cellfun(@(k)obj.elements{k}.name, num2cell(1:ne), ...
                            'UniformOutput', false );

            %% suppress output detail (default for large systems)
            suppress = mpopt.out.suppress_detail;
            if suppress == -1
                if obj.elements.is_index_name('bus') && ...
                        obj.elements.bus.nr > 500
                    suppress = 1;
                else
                    suppress = 0;
                end
            end

            %% element struct with all fields == 0, names are dm element names
            s0 = cell2struct(num2cell(zeros(1,ne)), names, 2);
            s1 = cell2struct(num2cell(ones( 1,ne)), names, 2);
            s2 = struct('all', 1, 'any', 1);

            out = struct( ...
                'all', mpopt.out.all, ...
                'any', 0, ...
                'sec', struct( ...
                    'cnt', s2, ...
                    'sum', s2, ...
                    'ext', s2, ...
                    'det', struct( ...
                        'all', -1, ...
                        'any', 0, ...
                        'elm', s1 ) ) );
            out.sec.cnt.all = out.all == 1 || (out.all == -1 && mpopt.out.sys_sum);
            out.sec.sum.all = out.sec.cnt.all;
            out.sec.ext.all = out.sec.cnt.all;
            out.sec.cnt.any = out.sec.cnt.all;
            out.sec.sum.any = out.sec.cnt.all;
            out.sec.ext.any = out.sec.cnt.all;
%             out.area_sum = out.all == 1 || (out.all == -1 && ~suppress && mpopt.out.area_sum);

            %% update detail options
            for k = 1:ne
                out.sec.det.elm.(names{k}) = out.all == 1 || ...
                    ( out.all == -1 && ~suppress && ...
                        (~isfield(mpopt.out, names{k}) || mpopt.out.(names{k})) );
            end
            out.sec.det.any = any(cell2mat( ...
                    cellfun(@double, struct2cell(out.sec.det.elm), ...
                            'UniformOutput', false) ));

            %% update any field
            out.any = out.all == 1 || ...
                        (out.all == -1 && (out.sec.cnt.all || out.sec.det.any ));

            %% return additional data
            if nargout > 1
                add = struct('s1', s1, 's0', s0, 'suppress', suppress, ...
                    'names', names, 'ne', ne);
            end
        end

        function h = pp_section_label(obj, label, blank_line)
            if nargin < 3
                blank_line = 1; %% include a blank line before the section label
            end

            width = 80;
            line1 = repmat('=', 1, width);
            line2 = ['|     ' label repmat(' ', 1, width-length(label)-7) '|'];
            if blank_line
                h = {'', line1, line2, line1};
            else
                h = {line1, line2, line1};
            end
        end

        function [obj, out_] = pretty_print(obj, mpopt, fd)
            if nargin < 3
                fd = 1;     %% print to stdio by default
                if nargin < 2
                    mpopt = mpoption();
                end
            end

            %% get output flags
            out = obj.pp_flags(mpopt);

            if out.any
                sections = obj.pp_section_list(out);  %% e.g. cnt, sum, ext, det, etc.
                for s = 1:length(sections)
                    out_s = out.sec.(sections{s});
                    if out_s.any && obj.pp_have_section(sections{s}, mpopt)
                        obj.pp_section(sections{s}, out_s, mpopt, fd);
                    end
                end
            end

            %% return output control flags
            if nargout > 1
                out_ = out;
            end
        end

        function sections = pp_section_list(obj, out)
            sections = {'cnt', 'sum', 'ext', 'det'};
        end

        function TorF = pp_have_section(obj, section, mpopt)
            %% section per-element info
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                TorF = dme.pp_have_section(section, mpopt, struct());
                if TorF
                    break;
                end
            end
        end

        function obj = pp_section(obj, section, out_s, mpopt, fd)
            %% section title & headers
            h = obj.pp_get_headers(section, out_s, mpopt);
            for k = 1:length(h)
                fprintf(fd, '%s\n', h{k});
            end

            %% section system info
            obj.pp_data(section, out_s, mpopt, fd);

            %% section per-element info
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                if out_s.all == -1
                    out_e = out_s.elm.(dme.name);
                else
                    out_e = out_s.all;
                end
                pp_args = struct('model', mpopt.model);
                dme.pretty_print(obj, section, out_e, mpopt, fd, pp_args);
            end
        end

        function h = pp_get_headers(obj, section, out_s, mpopt)
            switch section
                case 'cnt'
                    h = obj.pp_get_headers_cnt(out_s, mpopt);
                case 'sum'
                    h = {''};
                case 'ext'
                    h = obj.pp_get_headers_ext(out_s, mpopt);
                case 'det'
                    h = {};
                otherwise
                    h = obj.pp_get_headers_other(section, out_s, mpopt);
            end
        end

        function h = pp_get_headers_cnt(obj, out_s, mpopt)
            h = [ obj.pp_section_label('System Summary', 0) ...
                 {  '  elements                on     off    total', ...
                    ' --------------------- ------- ------- -------' } ];
        end

        function h = pp_get_headers_ext(obj, out_s, mpopt)
            h = {   '', ...
                    '                                           minimum                        maximum', ...
                    '                               -----------------------------  -----------------------------' };
        end

        function obj = pp_data(obj, section, out_s, mpopt, fd)
        end

        function obj = set_bus_v_lims_via_vg(obj, use_vg)
            bus_dme = obj.elements.bus;
            gen_dme = obj.elements.gen;
            gbus = bus_dme.i2on(gen_dme.bus(gen_dme.on));   %% buses of online gens
            nb = bus_dme.n;
            ng = gen_dme.n;

            %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
            Cg = sparse(gbus, (1:ng)', 1, nb, ng);
            Vbg = Cg * sparse(1:ng, 1:ng, gen_dme.vm_setpoint, ng, ng);
            vm_ub = max(Vbg, [], 2);    %% zero for non-gen buses, else max vm_setpoint of gens @ bus
            ib = find(vm_ub);               %% buses with online gens
            vm_lb = max(2*Cg - Vbg, [], 2); %% same as vm_ub, except min vm_setpoint of gens @ bus
            vm_lb(ib) = 2 - vm_lb(ib);

            if use_vg == 1      %% use vm_setpoint setpoint directly
                bus_dme.vm_ub(ib) = vm_ub(ib);  %% ub set by max vm_setpoint @ bus
                bus_dme.vm_lb(ib) = vm_lb(ib);  %% lb set by min vm_setpoint @ bus
                bus_dme.vm_start(ib) = vm_ub(ib);
            elseif use_vg > 0 && use_vg < 1     %% fractional value
                %% use weighted avg between original vm_lb/vm_ub limits and vm_setpoint
                bus_dme.vm_ub(ib) = (1-use_vg) * bus_dme.vm_ub(ib) + use_vg * vm_ub(ib);
                bus_dme.vm_lb(ib) = (1-use_vg) * bus_dme.vm_lb(ib) + use_vg * vm_lb(ib);
            else
                error('mp.data_model/set_bus_v_lims_via_vg: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
            end

            %% update bus table as well (for output)
            bus_dme.tab.vm_ub(bus_dme.on(ib)) = bus_dme.vm_ub(ib);
            bus_dme.tab.vm_lb(bus_dme.on(ib)) = bus_dme.vm_lb(ib);
        end
    end     %% methods
end         %% classdef
