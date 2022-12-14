classdef (Abstract) nm_element < handle
%MP.NM_ELEMENT  Abstract base class for MATPOWER network model elements
%   NME = MP.NM_ELEMENT()
%
%   Each concrete subclass must also inherit from a subclass of MP.FORM.
%
%   Properties
%       np : number of ports per element
%       nn : number of nodes per element (created by this element type)
%       nz : number of non-voltage state variables per element
%       nk : number of elements
%       C : cell array of sparse element-node incidence matrices,
%           where C{j}(i,k) is 1 if port j of element k is connected to node i
%       D : cell array of sparse incidence matrices for Z variables,
%           where D{j}(i,k) is 1 if j-th Z variable for element k corresponds
%           to element i of system Z
%
%   Methods
%       name() - name of element type (constant across formations)
%       count() - returns the number of elements of this type in dm, sets nme.nk
%       get_nv_()
%       x2vz()
%       incidence_matrix()
%       display()
%
%   Abstract Methods
%       add_nodes()
%       add_states()
%       add_vvars()
%       add_zvars()
%       build_params() - build model parameters from data model

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nk = 0;     %% number of elements of this type loaded
        C = [];     %% stacked element-node incidence matrices,
                    %% where C(i,kk) is 1 if port j of element k is
                    %% connected to node i, and kk = k + (j-1)*np
        D = [];     %% stacked sparse incidence matrices for
                    %% Z variables, where D(i,kk) is 1 if z-variable j
                    %% of element k is the i-th system z-variable
                    %% and kk = k + (j-1)*nz
        soln        %% struct for storing solved states, quantities
    end

    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function np = np(obj)
            np = 0;     %% number of ports per element
        end

        function nn = nn(obj)
            nn = 0;     %% number of nodes per element (created by element)
        end

        function nz = nz(obj)
            nz = 0;     %% number of (possibly complex) non-voltage states per element
        end

        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function mme = math_model_element(obj, mm, name)
            if nargin < 3
                name = obj.name;
            end
            if mm.elements.is_index_name(name)
                mme = mm.elements.(name);
            else
                mme = [];
            end
        end

        function nk = count(obj, dm)
            nk = dm.online(obj.name);
            obj.nk = nk;    %% update the count stored internally
        end

        function obj = add_nodes(obj, nm, dm)
            if obj.nn == 1
                nm.add_node(obj.name, obj.nk);
            elseif obj.nn > 1
                nm.init_indexed_name('node', obj.name, {obj.nn});
                for k = 1:obj.nn
                    nm.add_node(obj.name, {k}, obj.nk);
                end
            end
%             nm.add_node(obj.name, obj.nn * obj.nk);
        end

        function obj = add_states(obj, nm, dm)
            if obj.nz == 1
                nm.add_state(obj.name, obj.nk);
            elseif obj.nz > 1
                nm.init_indexed_name('state', obj.name, {obj.nz});
                for k = 1:obj.nz
                    nm.add_state(obj.name, {k}, obj.nk);
                end
            end
%             nm.add_state(obj.name, obj.nz * obj.nk);
        end

        function obj = add_vvars(obj, nm, dm, idx)
        end

        function obj = add_zvars(obj, nm, dm, idx)
        end

        function obj = build_params(obj, nm, dm)
            %% construct incidence matrices
            if obj.np
                nidx = obj.node_indices(nm, dm);
                obj.C = obj.incidence_matrix(nm.getN('node'), nidx{:});
            end
            if obj.nz
                sidx = nm.get_state_idx(obj.name);
                if iscell(sidx)
                    obj.D = obj.incidence_matrix(nm.getN('state'), sidx{:});
                else
                    obj.D = obj.incidence_matrix(nm.getN('state'), sidx);
                end
            else
                obj.D = obj.incidence_matrix(nm.getN('state'));
            end
        end

        function nv_ = get_nv_(obj, sysx);
            % sysx : 1 = system x_, 0 = element class x_

            %% get sizes
            if sysx
                nv_ = size(obj.C, 1);
%                 nz_ = size(obj.D, 1);
            else
                nv_ = obj.nk * obj.np;
%                 nz_ = obj.nk * obj.nz;
            end
        end

        function [v_, z_, vi_] = x2vz(obj, x_, sysx, idx);
            % sysx : 1 = system x_, 0 = element class x_
            % if x_ is a matrix, each output will have the same number of
            % columns, each column considered a separate instance of the vectors

            %% split x_
            nv_ = obj.get_nv_(sysx);
            v_ = x_(1:nv_, :);
            z_ = x_(nv_+1:end, :);

            %% set full port voltages and states for element class
            if sysx         %% system x_ is provided, convert to x_ for ports
                v_ = obj.C' * v_;   %% full port voltages for element class
                z_ = obj.D' * z_;   %% full states for element class
            end

            %% port voltages for selected ports
            if nargout > 2
                if isempty(idx)
                    vi_ = v_;
                else
                    vi_ = v_(idx, :);
                end
            end
        end

        function nidxs = node_indices(obj, nm, dm, cxn_type, cxn_idx_prop, cxn_type_prop)
            %% nidxs = obj.node_indices(nm, dm, dme, cxn_type, cxn_idx_prop)
            %% nidxs = obj.node_indices(nm, dm, dme, cxn_type, cxn_idx_prop, cxn_type_prop)
            %% A connection (cxn) is a mapping of a set of ports of element
            %% of type A (e.g. 'branch') to set of nodes created by elements
            %% of types B1, B2, etc. (e.g. 'bus'). We call the node creating
            %% elements "junction" (jxn) elements. A single connection links
            %% each type A element to exactly one type B element, where each
            %% link consists of N ports from A and N nodes from B, and N is
            %% determined by the number of nodes created by each type B element.
            %% Each of the following arguments can be a char array or cell
            %% array of char arrays.
            %%  cxn_type - name(s) of element type of junction elements to
            %%      which this element has connections, e.g. 'bus'
            %%      3 options:
            %%          1. single type for all connections, cxn_type is
            %%              a single char array, cxn_type_prop is empty
            %%          2. each connection has its own type, cxn_type is
            %%              a cell array of same dimension as cxn_idx_prop,
            %%              cxn_type_prop is empty
            %%          3. each individual element has it's own type,
            %%              cxn_type is a cell array and cxn_type_prop
            %%              provides the index into cxn_type for each element
            %%  cxn_idx_prop - name(s) of DME property(ies) containing
            %%      indices of the junction elements defining the connections,
            %%      e.g. {'bus_fr', 'bus_to'}; an empty char array signifies
            %%      a connection from the element to itself, i.e its ports
            %%      are connected to the nodes it created
            %%  cxn_type_prop - name(s) of DME property(ies) containing
            %%      type indices of the junction elements defining the
            %%      connections, where the type indices are indices into the
            %%      cxn_type cell array

            dme = obj.data_model_element(dm);   %% data model element for obj

            if nargin < 4
                cxn_type      = dme.cxn_type;
                cxn_idx_prop  = dme.cxn_idx_prop;
                cxn_type_prop = dme.cxn_type_prop;
            elseif nargin < 6
                cxn_type_prop = '';
            end

            if ischar(cxn_idx_prop)
                cxn_idx_prop = { cxn_idx_prop };
            end
            if ischar(cxn_type_prop) && ~isempty(cxn_type_prop)
                cxn_type_prop = { cxn_type_prop };
            end
            nc = length(cxn_idx_prop);  %% number of connections per element
            if iscell(cxn_type)
                if isempty(cxn_type_prop)
                    assert(length(cxn_type) == nc, 'mp.nm_element/node_indices: requires CXN_TYPE_PROP if specifying more than one CXN_TYPE and CXN_TYPE length does not match that of CXN_IDX_PROP');
                else
                    assert(length(cxn_type_prop) == nc, 'mp.nm_element/node_indices: CXN_TYPE_PROP length must match that of CXN_IDX_PROP');
                end
            end

            nidxs = cell(1, obj.np);    %% init node indices for each port

            if isempty(cxn_type_prop)
                if ischar(cxn_type)     %%-----  one type for all connections
                    jxn_dme = dm.elements.(cxn_type);   %% data model element for cxn_type
                    jxn_nme = nm.elements.(cxn_type);   %% corresp. net model element
                    ni = nm.get_node_idx(cxn_type);     %% corresponding node indices

                    %% number of connections
                    assert(obj.np == nc * jxn_nme.nn, ...
                        'mp.nm_element/node_indices: number of %s connections per %s (%d) x number of nodes per %s (%d) should equal the number of ports per %s (%d)', ...
                        cxn_type, obj.name, nc, cxn_type, jxn_nme.nn, obj.name, obj.np);

                    p = 1;                      %% initialize port index
                    for k = 1:nc                %% for connection k
                        %% jxn indices for cxn k of online elements
                        if isempty(cxn_idx_prop{k})
                            j = dme.on;
                        else
                            j = dme.(cxn_idx_prop{k})(dme.on);
                        end
                        jon = jxn_dme.i2on(j);  %% indices into online jxn element
                        if iscell(ni)
                            for n = 1:jxn_nme.nn
                                nidxs{p} = ni{n}(jon);
                                p = p + 1;          %% increment port index
                            end
                        else
                            nidxs{p} = ni(jon); %% node indices for single port of cxn k
                            p = p + 1;          %% increment port index
                        end
                    end
                else                    %%-----  one type for each connection
                    p = 1;                      %% initialize port index
                    for k = 1:nc                %% for connection k
                        jxn_dme = dm.elements.(cxn_type{k});%% data model element for cxn_type{k}
                        jxn_nme = nm.elements.(cxn_type{k});%% corresp. net model element
                        ni = nm.get_node_idx(cxn_type{k});  %% corresponding node indices

                        %% number of connections
                        assert(p <= obj.np, ...
                            'mp.nm_element/node_indices: port index for %s connections for %s (%d) exceeds the available number of ports per %s (%d)', ...
                            cxn_type{k}, obj.name, p, obj.name, obj.np);

                        %% jxn indices for cxn k of online elements
                        j = dme.(cxn_idx_prop{k})(dme.on);
                        jon = jxn_dme.i2on(j);  %% indices into online jxn element
                        if iscell(ni)
                            for n = 1:jxn_nme.nn
                                nidxs{p} = ni{n}(jon);
                                p = p + 1;          %% increment port index
                            end
                        else
                            nidxs{p} = ni(jon); %% node indices for single port of cxn k
                            p = p + 1;          %% increment port index
                        end
                    end
                end
            else                        %%-----  individual types per element/connection
                nidx = cell(1, nc);         %% init node indices for each cxn
                [nidx{:}] = deal(zeros(dme.n, 1));

                for i = 1:length(cxn_type)
                    if dm.elements.is_index_name(cxn_type{i})
                        jxn_dme = dm.elements.(cxn_type{i});    %% data model element for cxn_type{i}
                        jxn_nme = nm.elements.(cxn_type{i});    %% corresp. net model element
                        ni = nm.get_node_idx(cxn_type{i});      %% corresponding node indices

                        %% number of connections
                        assert(obj.np == nc * jxn_nme.nn, ...
                            'mp.nm_element/node_indices: number of %s connections per %s (%d) x number of nodes per %s (%d) should equal the number of ports per %s (%d)', ...
                            cxn_type{i}, obj.name, nc, cxn_type{i}, jxn_nme.nn, obj.name, obj.np);

                        for k = 1:nc                %% for connection k
                            %% for cxn k of online elements
                            j  = dme.(cxn_idx_prop{k})(dme.on); %% jxn indices
                            jt = dme.(cxn_type_prop{k})(dme.on);%% jxn types
                            %% assign node indices for 1st port of cxn k, jxn type i
                            nidx{k}(jt == i) = ni(jxn_dme.i2on(j(jt == i)));
                        end
                    end
                end
                p = 1;                      %% initialize port index
                for k = 1:nc                %% for connection k
                    for n = 0:jxn_nme.nn-1
                        nidxs{p} = nidx{k} + n * dme.n; %% node indices for port p
                        p = p + 1;          %% increment port index
                    end
                end
            end
        end

        function CD = incidence_matrix(obj, m, varargin)
            %% obj.incidence_matrix(m, idx1, idx2, ...)
            %% m = total number of nodes / states
            n = length(varargin);   %% number of ports/z-vars
            if n == 1
                CD = sparse(varargin{1}, 1:obj.nk, 1, m, obj.nk);
            elseif n > 1
                blocks = cell(1, n);
                for i = 1:n
                    blocks{i} = sparse(varargin{i}, 1:obj.nk, 1, m, obj.nk);
                end
                CD = horzcat(blocks{:});
            else
                CD = sparse(m, 0);
            end
        end

%         function A = getA(obj, tr)
%             if nargin < 1
%                 C = obj.C';
%                 D = obj.D';
%             else
%                 C = obj.C;
%                 D = obj.D;
%             end
%             [mC, nC] = size(C);
%             [mD, nD] = size(D);
%             A = [ C sparse(mC, nD); sparse(mD, nC) D ];
%         end
% 
%         function Ap = getAprime(obj, tr)
%             if nargin < 1
%                 C = obj.C';
%                 D = obj.D';
%             else
%                 C = obj.C;
%                 D = obj.D;
%             end
%             [mC, nC] = size(C);
%             [mD, nD] = size(D);
%             Ap = [   C sparse(mC, nC+2*nD);
%                     sparse(mC, nC) C sparse(mC, 2*nD);
%                     sparse(mD, 2*nC) D sparse(mD, nD);
%                     sparse(mD, 2*nC+nD) D ];
%         end

        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('NETWORK MODEL ELEMENT NAME  : %s\n', obj.name);
            fprintf('NETWORK MODEL ELEMENT CLASS : %s\n', class(obj));
            fprintf('    # OF ELEMENTS           : %d\n', obj.nk);
            fprintf('    # OF NODES/ELEM         : %d\n', obj.nn);
            fprintf('    # OF PORTS/ELEM         : %d\n', obj.np);
            fprintf('    # OF NON-V STATES/ELEM  : %d\n', obj.nz);
            if isa(obj, 'mp.form')
                fprintf('    FORMULATION NAME        : %s\n', obj.form_name());
                fprintf('    FORMULATION TAG         : %s\n', obj.form_tag());
                fprintf('    FORMULATION CLASS       : %s\n', obj.find_form_class());
                fprintf('    MODEL PARAMETERS');
                model_params = obj.model_params();
                for j = 1:length(model_params)
                    pn = model_params{j};   %% parameter name
                    if j == 1
                        fmt = '%7s : ';
                    else
                        fmt = '%27s : ';
                    end
                    if isempty(obj.(pn))
                        fprintf([fmt '-\n'], pn);
                    else
                        [m, n] = size(obj.(pn));
                        if ~full(any(any(obj.(pn))))
                            s = '(all zeros)';
                        else
                            s = '';
                        end
                        fprintf([fmt '%d x %-7d%s\n'], pn, m, n, s);
                    end
                end
            end
        end
    end     %% methods
end         %% classdef
