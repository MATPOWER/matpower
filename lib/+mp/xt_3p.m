classdef xt_3p < mp.extension
% mp.xt_3p - |MATPOWER| extension to add unbalanced three-phase elements.
%
% For AC power flow, continuation power flow, and optimial power flow problems,
% adds six new element types:
%
%   - ``'bus3p'`` - 3-phase bus
%   - ``'gen3p'`` - 3-phase generator
%   - ``'load3p'`` - 3-phase load
%   - ``'line3p'`` - 3-phase distribution line
%   - ``'xfmr3p'`` - 3-phase transformer
%   - ``'buslink'`` - 3-phase to single phase linking element
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
% mp.xt_3p Methods:
%   * dmc_element_classes - add six classes to data model converter elements
%   * dm_element_classes - add six classes to  data model elements
%   * nm_element_classes - add six classes to network model elements
%   * mm_element_classes - add six classes to mathematical model elements
%
% See the :ref:`sec_customizing` and :ref:`sec_extensions` sections in the
% |MATPOWER-Dev-Manual| for more information, and specifically the
% :ref:`sec_element_classes` section and the :ref:`tab_element_class_modifiers`
% table for details on *element class modifiers*.
%
% See also mp.extension.

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
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            % Add six classes to data model converter elements.
            %
            % For ``'mpc2`` data formats, adds the classes:
            %
            %   - mp.dmce_bus3p_mpc2
            %   - mp.dmce_gen3p_mpc2
            %   - mp.dmce_load3p_mpc2
            %   - mp.dmce_line3p_mpc2
            %   - mp.dmce_xfmr3p_mpc2
            %   - mp.dmce_buslink_mpc2

            switch fmt
                case 'mpc2'
                    dmc_elements = { ...
                        @mp.dmce_bus3p_mpc2, @mp.dmce_gen3p_mpc2, ...
                        @mp.dmce_load3p_mpc2, @mp.dmce_line3p_mpc2, ...
                        @mp.dmce_xfmr3p_mpc2, @mp.dmce_buslink_mpc2 ...
                    };
                otherwise
                    dmc_elements = {};
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            % Add six classes to data model elements.
            %
            % For ``'PF'`` and ``'CPF'`` tasks, adds the classes:
            %
            %   - mp.dme_bus3p
            %   - mp.dme_gen3p
            %   - mp.dme_load3p
            %   - mp.dme_line3p
            %   - mp.dme_xfmr3p
            %   - mp.dme_buslink
            %
            % For ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.dme_bus3p_opf
            %   - mp.dme_gen3p_opf
            %   - mp.dme_load3p_opf
            %   - mp.dme_line3p_opf
            %   - mp.dme_xfmr3p_opf
            %   - mp.dme_buslink_opf

            switch task_tag
                case {'PF', 'CPF'}
                    dm_elements = { ...
                        @mp.dme_bus3p, @mp.dme_gen3p, @mp.dme_load3p, ...
                        @mp.dme_line3p, @mp.dme_xfmr3p, @mp.dme_buslink ...
                    };
                case 'OPF'
                    dm_elements = { ...
                        @mp.dme_bus3p_opf, @mp.dme_gen3p_opf, ...
                        @mp.dme_load3p_opf, @mp.dme_line3p_opf, ...
                        @mp.dme_xfmr3p_opf, @mp.dme_buslink_opf ...
                    };
                otherwise
                    dm_elements = {};
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            % Add six classes to network model elements.
            %
            % For *cartesian* voltage formulations, adds the classes:
            %
            %   - mp.nme_bus3p_acc
            %   - mp.nme_gen3p_acc
            %   - mp.nme_load3p
            %   - mp.nme_line3p
            %   - mp.nme_xfmr3p
            %   - mp.nme_buslink_acc
            %
            % For *polar* voltage formulations, adds the classes:
            %
            %   - mp.nme_bus3p_acp
            %   - mp.nme_gen3p_acp
            %   - mp.nme_load3p
            %   - mp.nme_line3p
            %   - mp.nme_xfmr3p
            %   - mp.nme_buslink_acp

            switch task_tag
                case {'PF', 'CPF'}
                    v_cartesian = mpopt.pf.v_cartesian;
                case {'OPF'}
                    v_cartesian = mpopt.opf.v_cartesian;
            end
            switch upper(mpopt.model)
                case 'AC'
                    if v_cartesian
                        nm_elements = { ...
                            @mp.nme_bus3p_acc, @mp.nme_gen3p_acc, ...
                            @mp.nme_load3p, @mp.nme_line3p, ...
                            @mp.nme_xfmr3p, @mp.nme_buslink_acc ...
                        };
                    else
                        nm_elements = { ...
                            @mp.nme_bus3p_acp, @mp.nme_gen3p_acp, ...
                            @mp.nme_load3p, @mp.nme_line3p, ...
                            @mp.nme_xfmr3p, @mp.nme_buslink_acp ...
                        };
                    end
                case 'DC'
                    nm_elements = {};       %% no modifications
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            % Add five classes to mathematical model elements.
            %
            % For ``'PF'`` and ``'CPF'`` tasks, adds the classes:
            %
            %   - mp.mme_bus3p
            %   - mp.mme_gen3p
            %   - mp.mme_line3p
            %   - mp.mme_xfmr3p
            %   - mp.mme_buslink_pf_acc *(cartesian)* or
            %     mp.mme_buslink_pf_acp *(polar)*
            %
            % For ``'OPF'`` tasks, adds the classes:
            %
            %   - mp.mme_bus3p_opf_acc *(cartesian)* or
            %     mp.mme_bus3p_opf_acp *(polar)*
            %   - mp.mme_gen3p_opf
            %   - mp.mme_line3p_opf
            %   - mp.mme_xfmr3p_opf
            %   - mp.mme_buslink_opf_acc *(cartesian)* or
            %     mp.mme_buslink_opf_acp *(polar)*

            switch task_tag
                case {'PF', 'CPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.pf.v_cartesian
                                mm_elements = { ...
                                    @mp.mme_bus3p, @mp.mme_gen3p, ...
                                    @mp.mme_line3p, @mp.mme_xfmr3p, ...
                                    @mp.mme_buslink_pf_acc ...
                                };
                            else
                                mm_elements = { ...
                                    @mp.mme_bus3p, @mp.mme_gen3p, ...
                                    @mp.mme_line3p, @mp.mme_xfmr3p, ...
                                    @mp.mme_buslink_pf_acp ...
                                };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                    end
                case {'OPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.opf.v_cartesian
                                mm_elements = { ...
                                    @mp.mme_bus3p_opf_acc, ...
                                    @mp.mme_gen3p_opf, ...
                                    @mp.mme_line3p_opf, ...
                                    @mp.mme_xfmr3p_opf, ...
                                    @mp.mme_buslink_opf_acc ...
                                };
                            else
                                mm_elements = { ...
                                    @mp.mme_bus3p_opf_acp, ...
                                    @mp.mme_gen3p_opf, ...
                                    @mp.mme_line3p_opf, ...
                                    @mp.mme_xfmr3p_opf, ...
                                    @mp.mme_buslink_opf_acp ...
                                };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                    end
            end
        end
    end     %% methods
end         %% classdef
