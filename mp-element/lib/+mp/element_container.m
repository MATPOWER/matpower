classdef (Abstract) element_container < handle
%MP.ELEMENT_CONTAINER  Mix-in class to handle named/ordered element object array

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        element_classes %% cell array of function handles of
                        %% constructors for classes for individual
                        %% element types, filled by subclass constructor
        elements
    end     %% properties

    methods
        function obj = modify_element_classes(obj, class_list)
            %% each element in the class_list cell array is either:
            %%  1 - a handle to a constructor to be appended to
            %%      obj.element_classes, or
            %%  2 - a 2-element cell array {A,B} where A is a handle to
            %%      a constructor to replace any element E in the list for
            %%      which isa(E(), B) is true, i.e. B is a char array
            %%  3 - a char array B, where any element E in the list for
            %%      which isa(E(), B) is true is removed from the list
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
