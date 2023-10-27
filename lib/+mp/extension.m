classdef (Abstract) extension < handle
% mp.extension - Abstract base class for |MATPOWER| extensions.
%
% This class serves as the framework for the |*MATPOWER*| **extension** API,
% providing a way to bundle a set of class additions and modifications together
% into a single named package.
%
% By default the methods in this class do nothing, but they can be overridden
% to customize essentially any aspect of a |MATPOWER| run. The first 5 methods
% are used to modify the default classes used to construct the task, data
% model converter, data, network, and/or mathematical model objects. The last
% 4 methods are used to add to or modify the classes used to construct
% the elements for each of the container types.
%
% By convention, |MATPOWER| extension objects (or cell arrays of them) are
% named ``mpx`` and |MATPOWER| extension class names begin with ``mp.xt``.
%
% mp.extension Methods:
%   * task_class - return handle to constructor for task object
%   * dmc_class - return handle to constructor for data model converter object
%   * dm_class - return handle to constructor for data model object
%   * nm_class - return handle to constructor for network model object
%   * mm_class - return handle to constructor for mathematical object
%   * dmc_element_classes - return element class modifiers for data model converter elements
%   * dm_element_classes - return element class modifiers for data model elements
%   * nm_element_classes - return element class modifiers for network model elements
%   * mm_element_classes - return element class modifiers for mathematical model elements
%
% See the :ref:`sec_customizing` and :ref:`sec_extensions` sections in the
% |MATPOWER-Dev-Manual| for more information, and specifically the
% :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers`
% table for details on *element class modifiers*.
%
% Example |MATPOWER| extensions:
%
%   - mp.xt_reserves - adds fixed zonal reserves to OPF
%   - mp.xt_3p - adds example prototype unbalanced three-phase elements for
%     AC PF, CPF, and OPF
%
% See also mp.task, mp.dm_converter, mp.data_model, mp.net_model, mp.math_model,
% mp.dmc_element, mp.dm_element, mp.nm_element, mp.mm_element.

%   MATPOWER
%   Copyright (c) 2022-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function task_class = task_class(obj, task_class, mpopt)
            % Return handle to constructor for task object.
            % ::
            %
            %   task_class = mpx.task_class(task_class, mpopt)
            %
            % Inputs:
            %   task_class (function handle) : default task constructor
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   task_class (function handle) : updated task constructor
        end

        function dmc_class = dm_converter_class(obj, dmc_class, fmt, mpopt)
            % Return handle to constructor for data model converter object.
            % ::
            %
            %   dmc_class = mpx.dm_converter_class(dmc_class, fmt, mpopt)
            %
            % Inputs:
            %   dmc_class (function handle) : default data model converter constructor
            %   fmt (char array) : data format tag, e.g. ``'mpc2'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   dmc_class (function handle) : updated data model converter constructor
        end

        function dm_class = data_model_class(obj, dm_class, task_tag, mpopt)
            % Return handle to constructor for data model object.
            % ::
            %
            %   dm_class = mpx.data_model_class(dm_class, task_tag, mpopt)
            %
            % Inputs:
            %   dm_class (function handle) : default data model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   dm_class (function handle) : updated data model constructor
        end

        function nm_class = network_model_class(obj, nm_class, task_tag, mpopt)
            % Return handle to constructor for network model object.
            % ::
            %
            %   nm_class = mpx.network_model_class(nm_class, task_tag, mpopt)
            %
            % Inputs:
            %   nm_class (function handle) : default network model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   nm_class (function handle) : updated network model constructor
        end

        function mm_class = math_model_class(obj, mm_class, task_tag, mpopt)
            % Return handle to constructor for mathematical model object.
            % ::
            %
            %   mm_class = mpx.math_model_class(mm_class, task_tag, mpopt)
            %
            % Inputs:
            %   mm_class (function handle) : default math model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   mm_class (function handle) : updated math model constructor
        end

        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Return element class modifiers for data model converter elements.
            % ::
            %
            %   dmc_elements = mpx.dmc_element_classes(dmc_class, fmt, mpopt)
            %
            % Inputs:
            %   dmc_class (function handle) : data model converter constructor
            %   fmt (char array) : data format tag, e.g. ``'mpc2'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   dmc_elements (cell array) : element class modifiers (see
            %       :ref:`tab_element_class_modifiers` in the
            %       |MATPOWER-Dev-Manual|)

            dmc_elements = {};      %% no modifications
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            % Return element class modifiers for data model elements.
            % ::
            %
            %   dm_elements = mpx.dm_element_classes(dm_class, task_tag, mpopt)
            %
            % Inputs:
            %   dm_class (function handle) : data model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   dm_elements (cell array) : element class modifiers (see
            %       :ref:`tab_element_class_modifiers` in the
            %       |MATPOWER-Dev-Manual|)

            dm_elements = {};       %% no modifications
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            % Return element class modifiers for network model elements.
            % ::
            %
            %   nm_elements = mpx.nm_element_classes(nm_class, task_tag, mpopt)
            %
            % Inputs:
            %   nm_class (function handle) : network model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   nm_elements (cell array) : element class modifiers (see
            %       :ref:`tab_element_class_modifiers` in the
            %       |MATPOWER-Dev-Manual|)

            nm_elements = {};       %% no modifications
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Return element class modifiers for mathematical model elements.
            % ::
            %
            %   mm_elements = mpx.mm_element_classes(mm_class, task_tag, mpopt)
            %
            % Inputs:
            %   mm_class (function handle) : mathematical model constructor
            %   task_tag (char array) : task tag, e.g. ``'PF'``, ``'CPF'``, ``'OPF'``
            %   mpopt (struct) : |MATPOWER| options struct
            %
            % Output:
            %   mm_elements (cell array) : element class modifiers (see
            %       :ref:`tab_element_class_modifiers` in the
            %       |MATPOWER-Dev-Manual|)

            mm_elements = {};       %% no modifications
        end
    end     %% methods
end         %% classdef
