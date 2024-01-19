classdef dme_gen_node_test < mp.dme_gen_opf
%MP.DME_GEN_NODE_TEST  MATPOWER data model class for gen data for T_NODE_TEST

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nbet        %% number of bus element types
        bus_etv     %% bus element type vector (all gens), 1 = nld, 2 = ld
    end     %% properties

    methods
        function name = cxn_type(obj)
            name = {'bus_nld', 'bus_ld'};
        end

        function name = cxn_idx_prop(obj)
            name = 'bus';
        end

        function name = cxn_type_prop(obj)
            name = 'bus_etv';
        end

        function obj = initialize(obj, dm)
            initialize@mp.dm_element(obj, dm);  %% call parent

            %% get bus mapping info
            obj.nbet = length(obj.cxn_type);
            for k = obj.nbet:-1:1
                if dm.elements.has_name(obj.cxn_type{k})
                    bus_dme{k} = dm.elements.(obj.cxn_type{k});
                    b2i_k{k} = bus_dme{k}.ID2i;
                else
                    bus_dme{k} = [];
                    b2i_k{k} = [];
                end
                n(k) = length(b2i_k{k});
            end

            %% expand individual b2i mappings to be same dimension
            n_max = max(n);
            b2i = zeros(n_max, 1);
            for k = obj.nbet:-1:1
                if n(k) < n_max
                    b2i_k{k}(n_max, 1) = 0;
                end
                b2i = b2i + b2i_k{k};
            end

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
            obj.bus_etv = zeros(size(obj.bus));
            for k = 1:obj.nbet
                gk = find(b2i_k{k}(obj.tab.bus));
                if ~isempty(bus_dme{k})
                    obj.bus_etv(gk) = bus_dme{k}.bus_eti;
                end
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            for k = 1:obj.nbet
                if dm.elements.has_name(obj.cxn_type{k})
                    bus_dme = dm.elements.(obj.cxn_type{k});
                    bs = bus_dme.tab.status;    %% bus element status

                    %% update status of gens at isolated/offline buses
                    gk = find(obj.bus_etv == bus_dme.bus_eti);
                    obj.tab.status(gk) = obj.tab.status(gk) & bs(obj.bus(gk));
                end
            end

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);

            %% set bus_dme.vm_control for any online gens at PV buses
            for k = 1:obj.nbet
                if dm.elements.has_name(obj.cxn_type{k})
                    bus_dme = dm.elements.(obj.cxn_type{k});
                    gk = find(obj.bus_etv(obj.on) == bus_dme.bus_eti);
                    obj.bus_on(gk) = bus_dme.i2on(obj.bus(obj.on(gk)));
                    bt = bus_dme.type(obj.bus_on(gk));  %% bus type for online gens
                    i = find(bt == mp.NODE_TYPE.REF | bt == mp.NODE_TYPE.PV);
                    bus_dme.vm_control(obj.bus_on(gk(i))) = 1;
                end
            end
        end

        function obj = apply_vm_setpoint(obj, dm)
            %

            % set starting bus voltage, if bus is voltage-controlled
            for k = 1:obj.nbet
                if dm.elements.has_name(obj.cxn_type{k})
                    bus_dme = dm.elements.(obj.cxn_type{k});
                    gk = find(obj.bus_etv(obj.on) == bus_dme.bus_eti);
                    i = find(bus_dme.vm_control(obj.bus_on(gk)));
                    bus_dme.vm_start(obj.bus_on(gk(i))) = obj.vm_setpoint(gk(i));
                end
            end
        end
    end     %% methods
end         %% classdef
