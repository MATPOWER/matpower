classdef (Abstract) element_container < handle
% mp.element_container - Mix-in class to handle named/ordered element object array.
%
% Implements an element container that is used for |MATPOWER| model and
% data model converter objects. Provides the properties to store the
% constructors for each element and the elements themselves. Also provides
% a method to modify an existing set of element constructors.
%
% mp.element_container Properties:
%   * element_classes - cell array of element constructors
%   * elements - a mp.mapped_array to hold the element objects
%
% mp.element_container Methods:
%   * modify_element_classes - modify an existing set of element constructors
%
% See also mp.mapped_array.

%   MATPOWER
%   Copyright (c) 2020-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        % Cell array of function handles of constructors for individual
        % elements, filled by constructor of subclass.
        element_classes

        % A mapped array (mp.mapped_array) to hold the element objects
        % included inside this container object.
        elements
    end     %% properties

    methods
        function obj = modify_element_classes(obj, class_list)
            % Modify an existing set of element constructors.
            % ::
            %
            %   obj.modify_element_classes(class_list)
            %
            % Input:
            %   class_list (cell array) : list of **element class modifiers**,
            %       where each modifier is one of the following:
            %
            %         1.  a handle to a constructor to **append** to
            %             ``obj.element_classes``, *or*
            %         2.  a char array ``B``, indicating to **remove** any
            %             element ``E`` in the list for which ``isa(E(), B)``
            %             is ``true``, *or*
            %         3.  a 2-element cell array ``{A,B}`` where ``A`` is a
            %             handle to a constructor to **replace** any element ``E``
            %             in the list for which ``isa(E(), B)`` is ``true``,
            %             i.e. ``B`` is a char array
            %
            %       Also accepts a single element class modifier of type 1
            %       or 2 *(A single type 3 modifier has to be enclosed
            %       in a single-element cell array to keep it from being
            %       interpreted as a list of 2 modifiers)*.
            %
            % Can be used to modify the list of element constructors in the
            % :attr:`element_classes` property by appending, removing, or
            % replacing entries. See :numref:`tab_element_class_modifiers`
            % in the |MATPOWER-Dev-Manual| for more information.

            if ~iscell(class_list)
                class_list = {class_list};
            end
            ec = obj.element_classes;   %% list to be updated
            for k = 1:length(class_list)
                c = class_list{k};
                if iscell(c)        %% it's a 2-d cell array
                    i = find(cellfun(@(e)isa(e(), c{2}), ec));  %% find c{2}
                    if ~isempty(i)
                        ec{i} = c{1};                   %% replace with c{1}
                    end
                elseif ischar(c)
                    i = find(cellfun(@(e)isa(e(), c), ec));     %% find c
                    if ~isempty(i)
                        ec(i) = [];                     %% delete
                    end
                else                %% it's a single function handle
                    ec{end+1} = c;  %%      append it
                end
            end
            obj.element_classes = ec;
        end
    end     %% methods
end         %% classdef
