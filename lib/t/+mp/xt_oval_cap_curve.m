classdef xt_oval_cap_curve < mp.extension
% mp.xt_oval_cap_curve - |MATPOWER| extension for OPF with oval gen PQ capability curves.
%
% For OPF problems, this extension restricts the output of each generator
% to lie within the half-oval-shaped region centered at (PMIN, Q0) and
% passing though (PMAX, Q0), (PMIN, QMIN) and (PMIN, QMAX).
%
% mp.xt_oval_cap_curve Methods:
%   * mm_element_classes - replace a class in mathematical model elements
%
% See the :ref:`sec_customizing` and :ref:`sec_extensions` sections in the
% |MATPOWER-Dev-Manual| for more information, and specifically the
% :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers`
% table for details on *element class modifiers*.
%
% See also mp.extension, mp.mme_gen_opf_ac_oval.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Replace a class in mathematical model elements.
            %
            % For ``'OPF'`` tasks, replaces mp.gen_opf_ac with
            % mp.gen_opf_ac_oval.

            switch task_tag
                case {'OPF'}
                    mm_elements = { {@mp.mme_gen_opf_ac_oval, 'mp.mme_gen_opf_ac'} };
                otherwise
                    mm_elements = {};       %% no modifications
            end
        end
    end     %% methods
end         %% classdef
