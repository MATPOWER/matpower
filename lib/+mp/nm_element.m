classdef (Abstract) nm_element < handle
% mp.nm_element - Abstract base class for |MATPOWER| **network model element** objects.
%
% A network model element object encapsulates all of the network model
% parameters for a particular element type. All network model element classes
% inherit from mp.nm_element and also, like the container, from a
% formulation-specific subclass of mp.form. Each element type typically
% implements its own subclasses, which are further subclassed per formulation.
% A given network model element object contains the aggregate network model
% parameters for all online instances of that element type, stored in the set
% of matrices and vectors that correspond to the formulation.
%
% By convention, network model element variables are named ``nme`` and network
% model element class names begin with ``mp.nme``.
%
% mp.mm_element Properties:
%   * nk - number of elements of this type
%   * C - stacked sparse element-node incidence matrices
%   * D - stacked sparse incidence matrices for *z*-variables
%   * soln - struct for storing solved states, quantities
%
% mp.mm_element Methods:
%   * name - get name of element type, e.g. ``'bus'``, ``'gen'``
%   * np - number of ports per element of this type
%   * nn - number of nodes per element, created by this element type
%   * nz - number of non-voltage state variables per element of this type
%   * data_model_element - get the corresponding data model element
%   * math_model_element - get the corresponding math model element
%   * count - get number of online elements in ``dm``, set :attr:`nk`
%   * add_nodes - add nodes to network model
%   * add_states - add non-voltage states to network model
%   * add_vvars - add real-valued voltage variables to network object
%   * add_zvars - add real-valued non-voltage state variables to network object
%   * build_params - build model parameters from data model
%   * get_nv_ - get number of *(possibly complex)* voltage variables
%   * x2vz - get port voltages and non-voltage states from combined state vector
%   * node_indices - construct node indices from data model element connection info
%   * incidence_matrix - construct stacked incidence matrix from set of index vectors
%   * node_types - get node type information
%   * set_node_type_ref - make the specified node a reference node
%   * set_node_type_pv - make the specified node a PV node
%   * set_node_type_pq - make the specified node a PQ node
%   * display - display the network model element object
%
% See the :ref:`sec_nm_element` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.net_model.

%       build_params() - build model parameters from data model

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nk = 0;     % *(integer)* number of elements of this type

        % *(sparse integer matrix)* stacked element-node incidence matrices,
        % where ``C(i,kk)`` is 1 if port *j* of element *k* is connected to
        % node *i*, and ``kk = k + (j-1)*np``
        C = [];

        % *(sparse integer matrix)* stacked incidence matrices for
        % *z*-variables (non-voltage state variables), where ``D(i,kk)`` is 1
        % if *z*-variable *j* of element *k* is the *i*-th system *z*-variable
        % and ``kk = k + (j-1)*nz``
        D = [];
        soln        % *(struct)* for storing solved states, quantities
    end

    methods
        function name = name(obj)
            % Get name of element type, e.g. ``'bus'``, ``'gen'``.
            % ::
            %
            %   name = nme.name()
            %
            % Output:
            %   name (char array) : name of element type, must be a valid
            %       struct field name
            %
            % Implementation provided by an element type specific subclass.

            name = '';      %% e.g. 'bus', 'gen'
        end

        function np = np(obj)
            % Number of ports per element of this type.
            % ::
            %
            %   np = nme.np()
            %
            % Output:
            %   np (integer) : number of ports per element of this type

            np = 0;
        end

        function nn = nn(obj)
            % Number of nodes per element, created by this element type.
            % ::
            %
            %   nn = nme.nn()
            %
            % Output:
            %  nn (integer) : number of ports per element of this type

            nn = 0;
        end

        function nz = nz(obj)
            % Number of non-voltage state variables per element of this type.
            % ::
            %
            %   nz = nme.nz()
            %
            % Output:
            %   nz (integer) : number of non-voltage state variables per
            %       element of this type

            nz = 0;     %% number of (possibly complex) non-voltage states per element
        end

        function dme = data_model_element(obj, dm, name)
            % Get the corresponding data model element.
            % ::
            %
            %   dme = nme.data_model_element(dm)
            %   dme = nme.data_model_element(dm, name)
            %
            % Inputs:
            %   dm (mp.data_model) : data model object
            %   name (char array) : *(optional)* name of element type
            %       *(default is name of this object)*
            %
            % Output:
            %   dme (mp.dm_element) : data model element object

            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function mme = math_model_element(obj, mm, name)
            % Get the corresponding math model element.
            % ::
            %
            %   mme = nme.math_model_element(mm)
            %   mme = nme.math_model_element(mm, name)
            %
            % Inputs:
            %   mm (mp.math_model) : math model object
            %   name (char array) : *(optional)* name of element type
            %       *(default is name of this object)*
            %
            % Output:
            %   mme (mp.mm_element) : math model element object

            if nargin < 3
                name = obj.name;
            end
            if mm.elements.has_name(name)
                mme = mm.elements.(name);
            else
                mme = [];
            end
        end

        function nk = count(obj, dm)
            % Get number of online elements of this type in ``dm``, set :attr:`nk`.
            % ::
            %
            %   nk = nme.count(dm)
            %
            % Input:
            %   dm (mp.data_model) : data model object
            %
            % Output:
            %   nk (integer) : number of online elements of this type

            nk = dm.online(obj.name);
            obj.nk = nk;    %% update the count stored internally
        end

        function obj = add_nodes(obj, nm, dm)
            % Add nodes to network model for this element.
            % ::
            %
            %   nme.add_nodes(nm, dm)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Add nodes to the network model object, based on value *nn*
            % returned by nn(). Calls the network model's
            % :meth:`add_node() <mp.net_model.add_node>` *nn* times.

            if obj.nn == 1
                nm.add_node(obj.name, obj.nk);
            elseif obj.nn > 1
                nm.init_indexed_name('node', obj.name, {obj.nn});
                for k = 1:obj.nn
                    nm.add_node(obj.name, {k}, obj.nk);
                end
            end
        end

        function obj = add_states(obj, nm, dm)
            % Add non-voltage states to network model for this element.
            % ::
            %
            %   nme.add_states(nm, dm)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Add non-voltage states to the network model object, based on
            % value *nz* returned by nz(). Calls the network model's
            % :meth:`add_state() <mp.net_model.add_state>` *nz* times.

            if obj.nz == 1
                nm.add_state(obj.name, obj.nk);
            elseif obj.nz > 1
                nm.init_indexed_name('state', obj.name, {obj.nz});
                for k = 1:obj.nz
                    nm.add_state(obj.name, {k}, obj.nk);
                end
            end
        end

        function obj = add_vvars(obj, nm, dm, idx)
            % Add real-valued voltage variables to network object.
            % ::
            %
            %   nme.add_vvars(nm, dm, idx)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Add real-valued voltage variables (*v*-variables) to the network
            % model object, for each port. Implementation depends on the
            % specific formulation (i.e. subclass of mp.form).
            %
            % For example, consider an element with *np* ports and an AC
            % formulation with polar voltage representation. The actual port
            % voltages are complex, but this method would call the network
            % model's :meth:`add_var() <mp.net_model.add_var>` twice for each
            % port, once for the voltage angle variables and once for the
            % voltage magnitude variables.
            %
            % Implemented by a formulation-specific subclass.
        end

        function obj = add_zvars(obj, nm, dm, idx)
            % Add real-valued non-voltage state variables to network object.
            % ::
            %
            %   nme.add_zvars(nm, dm, idx)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   idx (cell array) : indices for named and indexed variables
            %
            % Add real-valued non-voltage state variables (*z*-variables) to
            % the network model object. Implementation depends on the specific
            % formulation (i.e. subclass of mp.form).
            %
            % For example, consider an element with *nz* z-variables and a
            % formulation in which these are complex. This method would call
            % the network model's :meth:`add_var() <mp.net_model.add_var>`
            % twice for each complex *z*-variable, once for the variables
            % representing the real part and once for the imaginary part.
            %
            % Implemented by a formulation-specific subclass.
        end

        function obj = build_params(obj, nm, dm)
            % Build model parameters from data model.
            % ::
            %
            %   nme.build_params(nm, dm)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %
            % Construction of incidence matrices :attr:`C` and :attr:`D` are
            % handled in this base class. Building of the formulation-specific
            % model parameters must be implemented by a formulation-specific
            % subclass. The subclass should call its parent in order to
            % construct the incidence matrices.
            %
            % See also incidence_matrix, node_indices.

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

        function nv_ = get_nv_(obj, sysx)
            % Get number of *(possibly complex)* voltage variables.
            % ::
            %
            %   nv_ = nme.get_nv_(sysx)
            %
            % Input:
            %   sysx (boolean) : if true the state ``x_`` refers to the full
            %       *(possibly complex)* system state *(all node voltages and
            %       system non-voltage states)*, otherwise it is the state
            %       vector for this specific element type *(port voltages and
            %       element non-voltage states)*
            %
            % Output:
            %   nv_ (integer) : number of *(possibly complex)* voltage variables
            %       in the state variable ``x_``, whose meaning depends on the
            %       ``sysx`` input

            %% get sizes
            if sysx
                nv_ = size(obj.C, 1);
%                 nz_ = size(obj.D, 1);
            else
                nv_ = obj.nk * obj.np;
%                 nz_ = obj.nk * obj.nz;
            end
        end

        function [v_, z_, vi_] = x2vz(obj, x_, sysx, idx)
            % Get port voltages and non-voltage states from combined state vector.
            % ::
            %
            %   [v_, z_, vi_] = nme.x2vz(x_, sysx, idx)
            %
            % Inputs:
            %   x_ (double) : *possibly complex* state vector
            %   sysx (boolean) : if true the state ``x_`` refers to the full
            %       *(possibly complex)* system state *(all node voltages and
            %       system non-voltage states)*, otherwise it is the state
            %       vector for this specific element type *(port voltages and
            %       element non-voltage states)*
            %   idx (integer) : vector of port indices of interest
            %
            % Outputs:
            %   v_ (double) : vector of *(possibly complex)* port voltages
            %   z_ (double) : vector of *(possibly complex)* non-voltage state
            %       variables
            %   vi_ (double) : vector of *(possibly complex)* port voltages for
            %       selected ports only, as indexed by ``idx``
            %
            % This method extracts voltage and non-voltage states from a
            % combined state vector, optionally with voltages for specific
            % ports only.
            %
            % Note, that this method can operate on multiple state vectors
            % simultaneously, by specifying ``x_`` as a matrix. In this case,
            % each output will have the same number of columns, one for each
            % column of the input ``x_``.

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
            % Construct node indices from data model element connection info.
            % ::
            %
            %   nidxs = nme.node_indices(nm, dm)
            %   nidxs = nme.node_indices(nm, dm, cxn_type, cxn_idx_prop)
            %   nidxs = nme.node_indices(nm, dm, cxn_type, cxn_idx_prop, cxn_type_prop)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   cxn_type (char array or cell array of char arrays) : name(s) of
            %       type(s) of junction elements, i.e. node-creating elements
            %       (e.g. ``'bus'``), to which this element connects; see
            %       mp.dm_element.cxn_type for more info
            %   cxn_idx_prop (char array or cell array of char arrays) : name(s)
            %       of property(ies) containing indices of junction elements
            %       that define connections (e.g. ``{'fbus', 'tbus'}``); see
            %       mp.dm_element.cxn_idx_prop for more info
            %   cxn_type_prop (char array or cell array of char arrays) :
            %       name(s) of properties containing type of junction elements
            %       for each connection, defaults to ``''`` if ``cxn_type`` and
            %       ``cxn_type_prop`` are provided, but not ``cxn_type_prop``;
            %       see mp.dm_element.cxn_type_prop for more info
            %
            % Output:
            %   nidxs (cell array) : 1 x *np* cell array of node index vectors
            %       for each port
            %
            % This method constructs the node index vectors for each port. That
            % is, element *p* of ``nidxs`` is the vector of indices of the
            % nodes to which port *p* of these elements are connected. These
            % node indices can be used to construct the element-node incidence
            % matrices that form :attr:`C`.
            %
            % By default, the connection information is obtained from the
            % corresponding data model element, as described in the
            % :ref:`sec_dm_element_cxn` section in the |MATPOWER-Dev-Manual|.
            %
            % See also incidence_matrix, mp.dm_element.cxn_type,
            % mp.dm_element.cxn_idx_prop, mp.dm_element.cxn_type_prop.

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
                    if dm.elements.has_name(cxn_type{i})
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
            % Construct stacked incidence matrix from set of index vectors.
            % ::
            %
            %   CD = nme.incidence_matrix(m, idx1, idx2, ...)
            %
            % Inputs:
            %   m (integer) : total number of nodes or states
            %   idx1 (integer) : index vector for nodes corresponding to
            %       this element's first port, or state variables corresponding
            %       to this element's first non-voltage state
            %   idx2 (integer) : same as ``idx1`` for second port or
            %       non-voltage state, and  so on
            %
            % Output:
            %   CD (sparse matrix) : stacked incidence matrix (``C`` for ports,
            %       ``D`` for states)
            %
            % Forms an *m* x *n* incidence matrix for each input index vector
            % ``idx``, where *n* is the dimension of ``idx``, and column ``j``
            % of the corresponding incidence matrix consists of all zeros with
            % a 1 in row ``idx(j)``.
            %
            % These incidence matrices are then stacked horizontally to form
            % a single matrix return value.

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

        function [ref, pv, pq] = node_types(obj, nm, dm, idx)
            % Get node type information.
            % ::
            %
            %   ntv           = nme.node_types(nm, dm)
            %   [ref, pv, pq] = nme.node_types(nm, dm)
            %             ... = nme.node_types(nm, dm, idx)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   idx (integer) : index *(not used in base method)*
            %
            % Outputs:
            %   ntv (integer) : node type vector, valid element values are:
            %
            %       - mp.NODE_TYPE.REF
            %       - mp.NODE_TYPE.PV
            %       - mp.NODE_TYPE.PQ
            %   ref (integer) : vector of indices of reference nodes
            %   pv (integer) : vector of indices of PV nodes
            %   pq (integer) : vector of indices of PQ nodes
            %
            % See also mp.NODE_TYPE.
            error('node_types() method not implemented for class ''%s''', class(obj));
        end

        function set_node_type_ref(obj, dm, idx)
            % Make the specified node a reference node.
            % ::
            %
            %   nme.set_node_type_ref(dm, idx)
            %
            % Inputs:
            %   dm (mp.data_model) : data model object
            %   idx (integer) : index of node to modify, this is the internal
            %       network model element index
            %
            % Set the specified node to type mp.NODE_TYPE.REF.
            %
            % Implementation provided by node-creating subclass.
            error('set_node_type_ref() method not implemented for class ''%s''', class(obj));
        end

        function set_node_type_pv(obj, dm, idx)
            % Make the specified node a PV node.
            % ::
            %
            %   nme.set_node_type_pv(dm, idx)
            %
            % Inputs:
            %   dm (mp.data_model) : data model object
            %   idx (integer) : index of node to modify, this is the internal
            %       network model element index
            %
            % Set the specified node to type mp.NODE_TYPE.PV.
            %
            % Implementation provided by node-creating subclass.
            error('set_node_type_pv() method not implemented for class ''%s''', class(obj));
        end

        function set_node_type_pq(obj, dm, idx)
            % Make the specified node a PQ node.
            % ::
            %
            %   nme.set_node_type_pq(dm, idx)
            %
            % Inputs:
            %   dm (mp.data_model) : data model object
            %   idx (integer) : index of node to modify, this is the internal
            %       network model element index
            %
            % Set the specified node to type mp.NODE_TYPE.PQ.
            %
            % Implementation provided by node-creating subclass.
            error('set_node_type_pq() method not implemented for class ''%s''', class(obj));
        end

        function display(obj)
            % Display the network model element object.
            %
            % This method is called automatically when omitting a semicolon
            % on a line that retuns an object of this class.
            %
            % Displays the details of the elements, including total number
            % of elements, nodes per element, ports per element, non-voltage
            % state per element, formulation name, tag, and class, and names
            % and dimensions of the model parameters.

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
