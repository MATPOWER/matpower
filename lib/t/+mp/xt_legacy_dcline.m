classdef xt_legacy_dcline < mp.extension
% mp.xt_legacy_dcline - |MATPOWER| extension to add legacy DC line elements.
%
% For AC and DC power flow, continuation power flow, and optimial power flow
% problems, adds a new element type:
%
%   - ``'legacy_dcline'`` - legacy DC line
%
% No changes are required for the task or container classes, so only the
% ``..._element_classes`` methods are overridden.
%
% The set of data model element classes depends on the task, with each OPF
% class inheriting from the corresponding class used for PF and CPF.
%
% The set of network model element classes depends on the formulation,
% specifically whether cartesian or polar representations are used for
% voltages.
%
% And the set of mathematical model element classes depends on both the task
% and the formulation.
%
% mp.xt_legacy_dcline Methods:
%   * dmc_element_classes - add a class to data model converter elements
%   * dm_element_classes - add a class to  data model elements
%   * nm_element_classes - add a class to network model elements
%   * mm_element_classes - add a class to mathematical model elements
%
% See the :ref:`sec_customizing` and :ref:`sec_extensions` sections in the
% |MATPOWER-Dev-Manual| for more information, and specifically the
% :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers`
% table for details on *element class modifiers*.
%
% See also mp.extension.

%   MATPOWER
%   Copyright (c) 2022-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Add a class to data model converter elements.
            %
            % For ``'mpc2`` data formats, adds the classes:
            %
            %   - mp.dmce_legacy_dcline_mpc2

            switch fmt
                case 'mpc2'
                    dmc_elements = { @mp.dmce_legacy_dcline_mpc2 };
                otherwise
                    dmc_elements = {};
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            % Add a class to data model elements.
            %
            % For ``'PF'`` and ``'CPF'`` tasks, adds the class:
            %
            %   - mp.dme_legacy_dcline
            %
            % For ``'OPF'`` tasks, adds the class:
            %
            %   - mp.dme_legacy_dcline_opf

            switch task_tag
                case {'PF', 'CPF'}
                    dm_elements = { @mp.dme_legacy_dcline };
                case 'OPF'
                    dm_elements = { @mp.dme_legacy_dcline_opf };
                otherwise
                    dm_elements = {};
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            % Add a class to network model elements.
            %
            % For DC formulations, adds the class:
            %
            %   - mp.nme_legacy_dcline_dc
            %
            % For AC *cartesian* voltage formulations, adds the class:
            %
            %   - mp.nme_legacy_dcline_acc
            %
            % For AC *polar* voltage formulations, adds the class:
            %
            %   - mp.nme_legacy_dcline_acp

            switch task_tag
                case {'PF', 'CPF'}
                    v_cartesian = mpopt.pf.v_cartesian;
                case {'OPF'}
                    v_cartesian = mpopt.opf.v_cartesian;
            end
            switch upper(mpopt.model)
                case 'AC'
                    if v_cartesian
                        nm_elements = { @mp.nme_legacy_dcline_acc };
                    else
                        nm_elements = { @mp.nme_legacy_dcline_acp };
                    end
                case 'DC'
                    nm_elements = { @mp.nme_legacy_dcline_dc };
                otherwise
                    nm_elements = {};
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Add a class to mathematical model elements.
            %
            % For ``'PF'`` and ``'CPF'`` tasks, adds the class:
            %
            %   - mp.mme_legacy_dcline_pf_dc *(DC formulation)* or
            %   - mp.mme_legacy_dcline_pf_ac *(AC formulation)*
            %
            % For ``'OPF'`` tasks, adds the class:
            %
            %   - mp.mme_legacy_dcline_opf_dc *(DC formulation)* or
            %   - mp.mme_legacy_dcline_opf_ac *(AC formulation)*

            switch task_tag
                case {'PF', 'CPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            mm_elements = { @mp.mme_legacy_dcline_pf_ac };
                        case 'DC'
                            mm_elements = { @mp.mme_legacy_dcline_pf_dc };
                    end
                case {'OPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            mm_elements = { @mp.mme_legacy_dcline_opf_ac };
                        case 'DC'
                            mm_elements = { @mp.mme_legacy_dcline_opf_dc };
                    end
                otherwise
                    dm_elements = {};
            end
        end
    end     %% methods
end         %% classdef
