classdef (Abstract) form < handle
% mp.form - Abstract base class for |MATPOWER| **formulation**.
%
% Used as a mix-in class for all **network model element** classes.
% That is, each concrete network model element class must inherit, at least
% indirectly, from both mp.nm_element and mp.form.
%
% mp.form provides properties and methods that are specific to the network
% model formulation (e.g. DC version, AC polar power version, etc.).
%
% For more details, see the :ref:`sec_net_model_formulations` section in the
% |MATPOWER-Dev-Manual| and the derivations in |TN5|.
% 
% mp.form Properties:
%    *subclasses provide properties for model parameters*
%
% mp.form Methods:
%   * form_name - get char array w/name of formulation
%   * form_tag - get char array w/short label of formulation
%   * model_params - get cell array of names of model parameters
%   * model_vvars - get cell array of names of voltage state variables
%   * model_zvars - get cell array of names of non-voltage state variables
%   * get_params - get network model element parameters
%   * find_form_class - get name of network element object's formulation subclass
%
% See also mp.nm_element.

%   MATPOWER
%   Copyright (c) 2019-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         mp_form_field = '';
%     end

    methods
        function name = form_name(obj)
            % Get user-readable name of formulation, e.g. ``'DC'``, ``'AC-cartesian'``, ``'AC-polar'``.
            % ::
            %
            %   name  = nme.form_name()
            %
            % Output:
            %   name (char array) : name of formulation
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.form.form_name: must be implemented in subclass');
        end

        function tag = form_tag(obj)
            % Get short label of formulation, e.g. ``'dc'``, ``'acc'``, ``'acp'``.
            % ::
            %
            %   tag  = nme.form_tag()
            %
            % Output:
            %   tag (char array) : short label of formulation
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.form.form_tag: must be implemented in subclass');
        end

        function params = model_params(obj)
            % Get cell array of names of model parameters.
            % ::
            %
            %   params  = nme.model_params()
            %
            % Output:
            %   params (cell array of char arrays) : names of object properies
            %       for model parameters
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.form.model_params: must be implemented in subclass');
        end

        function vtypes = model_vvars(obj)
            % Get cell array of names of voltage state variables.
            % ::
            %
            %   vtypes = nme.model_vvars()
            %
            % Output:
            %   vtypes (cell array of char arrays) : names of network object
            %       properties for voltage state variables
            %
            % The network model object, which inherits from mp_idx_manager,
            % uses these values as set types for tracking its voltage
            % state variables.
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.form.model_vvars: must be implemented in subclass');
        end

        function vtypes = model_zvars(obj)
            % Get cell array of names of non-voltage state variables.
            % ::
            %
            %   vtypes = nme.model_zvars()
            %
            % Output:
            %   vtypes (cell array of char arrays) : names of network object
            %       properties for voltage state variables
            %
            % The network model object, which inherits from mp_idx_manager,
            % uses these values as set types for tracking its non-voltage
            % state variables.
            %
            % *Note: This is an abstract method that must be implemented
            % by a subclass.*

            error('mp.form.model_zvars: must be implemented in subclass');
        end

        function varargout = get_params(obj, idx, names)
            % Get network model element parameters.
            % ::
            %
            %   [p1, p2, ..., pN] = nme.get_params(idx)
            %   pA = nme.get_params(idx, nameA)
            %   [pA, pB, ...] = nme.get_params(idx, {nameA, nameB, ...})
            %
            % Inputs:
            %   idx (integer) : *(optional)* vector of indices of ports of
            %       interest, if empty or missing, returns parameters
            %       corresponding to all ports
            %   names (char array or cell array of char arrays) : *(optional)*
            %       name(s) of parameters to return
            %
            % Outputs:
            %   p1, p2, ..., pN : full set of parameters in canonical order
            %   pA, pB : parameters specified by ``names``
            %
            % If a particular parameter in the object is empty, this method
            % returns a sparse zero matrix or vector of the appropriate size.

            if nargin < 3
                names = obj.model_params();
            end
            if nargin < 2
                idx = [];
            end
            if ~iscell(names)
                names = {names};
            end
            np = obj.nk * obj.np;
            ncols = [1; np; obj.nk * obj.nz];   %% 1 = 1, 2 = np, 3 = nz

            if nargout > length(names)
                namestr = sprintf(' ''%s''', names{:})
                error('mp.form.get_params: return values must correspond to%s', namestr);
            end

            varargout = cell(1, nargout);
            if isempty(idx)     %% all ports
                for k = 1:nargout
                    P = obj.(names{k});
                    if isempty(P)
                        varargout{k} = sparse(np, ncols(obj.param_ncols.(names{k})));
                    else
                        varargout{k} = P;
                    end
                end
            else                %% selected ports
                ni = length(idx);
                for k = 1:nargout
                    P = obj.(names{k});
                    if isempty(P)
                        varargout{k} = sparse(ni, ncols(obj.param_ncols.(names{k})));
                    else
                        varargout{k} = P(idx, :);
                    end
                end
            end
        end

        function form_class = find_form_class(obj)
            % Get name of network element object's formulation subclass.
            % ::
            %
            %   form_class = nme.find_form_class()
            %
            % Output:
            %   form_class (char array) - name of the mp.form subclass
            %
            % Selects from this netork model elements parent classes, the
            % mp.form subclass, that is not a subclass of mp.nm_element, with
            % the longest inheritance path back to mp.form.

            if isa(obj, 'mp.form')
                tab = obj.superclass_tab({'mp.form', 'mp.nm_element'});

                %% select among classes that are not mp.nm_element's ...
                j = find(tab.ii(:, 2) == 0);
                %% ... the mp.form class with longest inheritance path
                %% back to mp.form
                [~, k] = max(tab.ii(j, 1));

                assert(~isempty(k));
                form_class = tab.name{j(k)};
            else
                form_class = '<none>';
            end
        end
    end     %% methods

    methods (Access=protected)
        function [tab, ii] = superclass_tab(obj, roots, mcls, tab, ii, level)
            % Called recursively to build the inheritance tree.
            % ::
            %
            %   [tab, ii] = nme.superclass_tab(roots, mcls, tab, ii, level)
            %
            % Inputs:
            %   roots (cell array) : names of root classes to search for
            %   mcls (object) : meta class object (undocumented MATLAB and
            %       Octave) for the class of interest
            %   tab (struct) : with fields:
            %
            %       - ``name`` -- cell array of class names
            %       - ``ii`` -- matrix where each row corresponds to entry in name
            %         field, and column *j* corresponds to *j*-th entry in
            %         ``roots``, value is generations of inheritance
            %         required to reach that root
            %   ii (integer) : row of ``tab.ii`` corresponding to ``mcls`` on
            %       input *(before traversing parent classes)*
            %   level (integer) : indent level for this class when printing
            %         (i.e. when hard-coded ``'verbose'`` variable not 0)
            %
            % Outputs:
            %   tab (struct) : *see input*
            %   ii (struct) : *see input*
            %
            % Called by find_form_class().

            verbose = 0;
            if nargin < 6
                level = 0;
                if nargin < 5
                    ii = [];
                    if nargin < 4
                        tab = struct('name', {{}}, 'ii', []);
                        if nargin < 3
                            mcls = meta.class.fromName(class(obj));
                            if nargin < 2
                                roots = {};
                            end
                        end
                    end
                end
            end
            n = length(roots);
            if isempty(ii)
                ii = zeros(1, n);
            end
            done = zeros(1, n);
            i0 = ii;
            for k = 1:n
                if strcmp(mcls.Name, roots{k})
                    done(k) = 1;
                end
            end
            if verbose
                prefix = repmat(' ', 1, 4*level);
                fprintf('%s %s\n', prefix, mcls.Name);
            end
            if have_feature('octave')
                sclist = mcls.SuperClassList;
            else
                sclist = {};
                for k = 1:length(mcls.SuperclassList)
                    sclist{end+1} = mcls.SuperclassList(k);
                end
            end
            for k = 1:length(sclist)
                [tab, iii] = obj.superclass_tab(roots, sclist{k}, tab, i0, level+1);
                for k = 1:n
                    if done(k) == 1
                        ii(k) = 1;
                    elseif iii(k) && i0(k) < iii(k) + 1
                        ii(k) = iii(k) + 1;
                    end
                end
            end
            if verbose
                fmt = repmat(' %d', 1, n);
                fprintf(['%s' fmt '\n'], prefix, ii);
            end
            tab.name{end+1, 1} = mcls.Name;
            tab.ii(end+1, 1:n) = ii;
        end
    end     %% methods (Access=protected)
end         %% classdef
