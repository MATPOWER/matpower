classdef (Abstract) mm_element < handle
% mp.mm_element - Abstract base class for |MATPOWER| **mathematical model element** objects.
%
% A math model element object typically does not contain any data, but only
% the methods that are used to build the math model and update the
% corresponding data model element once the math model has been solved.
%
% All math model element classes inherit from mp.mm_element. Each element
% type typically implements its own subclasses, which are further subclassed
% where necessary per task and formulation, as with the container class.
%
% By convention, math model element variables are named ``mme`` and math model
% element class names begin with ``mp.mme``.
%
% mp.mm_element Methods:
%   * name - get name of element type, e.g. ``'bus'``, ``'gen'``
%   * data_model_element - get corresponding data model element
%   * network_model_element - get corresponding network model element
%   * add_vars - add math model variables for this element
%   * add_constraints - add math model constraints for this element
%   * add_costs - add math model costs for this element
%   * data_model_update - update the corresponding data model element
%   * data_model_update_off - update offline elements in corresponding data model element
%   * data_model_update_on - update online elements in corresponding data model element
%
% See the :ref:`sec_mm_element` section in the |MATPOWER-Dev-Manual| for more
% information.
%
% See also mp.math_model.

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
        function name = name(obj)
            % Get name of element type, e.g. ``'bus'``, ``'gen'``.
            % ::
            %
            %   name = mme.name()
            %
            % Output:
            %   name (char array) : name of element type, must be a valid
            %       struct field name
            %
            % Implementation provided by an element type specific subclass.

            name = '';
        end

        function dme = data_model_element(obj, dm, name)
            % Get corresponding data model element.
            % ::
            %
            %   dme = mme.data_model_element(dm)
            %   dme = mme.data_model_element(dm, name)
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

        function nme = network_model_element(obj, nm, name)
            % Get corresponding network model element.
            % ::
            %
            %   nme = mme.network_model_element(nm)
            %   nme = mme.network_model_element(nm, name)
            %
            % Inputs:
            %   nm (mp.net_model) : network model object
            %   name (char array) : *(optional)* name of element type
            %       *(default is name of this object)*
            %
            % Output:
            %   nme (mp.nm_element) : network model element object

            if nargin < 3
                name = obj.name;
            end
            nme = nm.elements.(name);
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
            % Add math model variables for this element.
            % ::
            %
            %   mme.add_vars(mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Implementation provided by a subclass.
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            % Add math model constraints for this element.
            % ::
            %
            %   mme.add_constraints(obj, mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Implementation provided by a subclass.
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            % Add math model costs for this element.
            % ::
            %
            %   mme.add_costs(obj, mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Implementation provided by a subclass.
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            % Update the corresponding data model element.
            % ::
            %
            %   mme.data_model_update(mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Call data_model_update_off() then data_model_update_on()
            % to update the data model for this element based on the
            % math model solution.
            %
            % See also data_model_update_off, data_model_update_on.

            %% update offline elements
            obj.data_model_update_off(mm, nm, dm, mpopt);

            %% update online elements
            obj.data_model_update_on(mm, nm, dm, mpopt);
        end

        function obj = data_model_update_off(obj, mm, nm, dm, mpopt)
            % Update offline elements in the corresponding data model element.
            % ::
            %
            %   mme.data_model_update_off(mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Set export variables for offline elements based on specs
            % returned by mp.dm_element.export_vars_offline_val.
            %
            % See also data_model_update, data_model_update_on.

            dme = obj.data_model_element(dm);
            if ~isempty(dme.off)
                vars = dme.export_vars();
                if ~isempty(vars)
                    s = dme.export_vars_offline_val();
                    for k = 1:length(vars)
                        if isfield(s, vars{k})
                            dme.tab.(vars{k})(dme.off) = s.(vars{k});
                        end
                    end
                end
            end
        end

        function obj = data_model_update_on(obj, mm, nm, dm, mpopt)
            % Update online elements in the corresponding data model element.
            % ::
            %
            %   mme.data_model_update_on(mm, nm, dm, mpopt)
            %
            % Inputs:
            %   mm (mp.math_model) : mathmatical model object
            %   nm (mp.net_model) : network model object
            %   dm (mp.data_model) : data model object
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Extract the math model solution relevant to this particular
            % element and update the corresponding data model element
            % for online elements accordingly.
            %
            % Implementation provided by a subclass.
            %
            % See also data_model_update, data_model_update_off.
        end
    end     %% methods
end         %% classdef
