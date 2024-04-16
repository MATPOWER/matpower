classdef dm_converter_mpc2_legacy < mp.dm_converter_mpc2
% mp.dm_converter_mpc2_legacy - Legacy |MATPOWER| **data model converter** for MATPOWER case v2.
%
% Adds to mp.dm_converter_mpc2 the ability to handle legacy user
% customization.
%
% mp.dm_converter_mpc2_legacy Methods:
%   * legacy_user_mod_inputs - pre-process legacy inputs for use-defined customization
%   * legacy_user_nln_constraints - pre-process legacy inputs for user-defined nonlinear constraints
%
% See also mp.dm_converter, mp.dm_converter_mpc2, mp.task_opf_legacy.

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
        function dm = legacy_user_mod_inputs(obj, dm, mpopt, dc)
            % Handle pre-processing of inputs related to legacy user-defined
            % variables, costs, and constraints. This includes optional
            % ``mpc`` fields ``A``, ``l``, ``u``, ``N``, ``fparm``, ``H1``,
            % ``Cw``, ``z0``, ``zl``, ``zu`` and ``user_constraints``.

            %% create (read-only) copies of individual fields for convenience
            mpc = dm.source;
            [~, bus, gen, ~, ~, A, l, u, mpopt, ...
                N, fparm, H, Cw, z0, zl, zu, ~] = opf_args(mpc, mpopt);

            %% data dimensions
            nb   = size(bus, 1);    %% number of buses
            ng   = size(gen, 1);    %% number of dispatchable injections
            nlin = size(A, 1);      %% number of linear user constraints
            nw = size(N, 1);        %% number of general cost vars, w
            if dc
                %% reduce A and/or N from AC dimensions to DC dimensions, if needed
                if nlin || nw
                    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
                    if nlin && size(A, 2) >= 2*nb + 2*ng
                        %% make sure there aren't any constraints on Vm or Qg
                        if any(any(A(:, acc)))
                            error('mp.dm_converter_mpc2_legacy.legacy_user_mod_inputs: attempting to solve DC OPF with user constraints on Vm or Qg');
                        end
                        A(:, acc) = [];         %% delete Vm and Qg columns
                    end
                    if nw && size(N, 2) >= 2*nb + 2*ng
                        %% make sure there aren't any costs on Vm or Qg
                        if any(any(N(:, acc)))
                            [ii, jj] = find(N(:, acc));
                            ii = unique(ii);    %% indices of w with potential non-zero cost terms from Vm or Qg
                            if any(Cw(ii)) || (~isempty(H) && any(any(H(:, ii))))
                                error('mp.dm_converter_mpc2_legacy.legacy_user_mod_inputs: attempting to solve DC OPF with user costs on Vm or Qg');
                            end
                        end
                        N(:, acc) = [];         %% delete Vm and Qg columns
                    end
                end
                nv = 0;
                nq = 0;
            else
                nv = nb;
                nq = ng;
            end

            %% get number of user vars, check consistency
            nx = nb+nv + ng+nq;  %% number of standard OPF control variables
            if nlin
                nz = size(A, 2) - nx;   %% number of user z variables
                if nz < 0
                    error('mp.dm_converter_mpc2_legacy.legacy_user_mod_inputs: user supplied A matrix must have at least %d columns.', nx);
                end
            else
                nz = 0;               %% number of user z variables
                if nw                 %% still need to check number of columns of N
                    if size(mpc.N, 2) ~= nx;
                        error('mp.dm_converter_mpc2_legacy.legacy_user_mod_inputs: user supplied N matrix must have %d columns.', nx);
                    end
                end
            end

            %% package up parameters
            z = struct( 'nz', nz, ...           %% num user variables
                        'z0', z0, ...
                        'zl', zl, ...
                        'zu', zu    );
            lin = struct(   'nlin', nlin, ...   %% num user linear constraints
                            'A', A, ...
                            'l', l, ...
                            'u', u  );
            cost =  struct( 'nw', nw, ...       %% num user cost rows
                            'N', N, ...
                            'Cw', Cw    );
            nlc = obj.legacy_user_nln_constraints(dm, mpopt);

            if ~isempty(H)
                cost.H = H;
            end
            if ~isempty(fparm)
                cost.dd = fparm(:, 1);
                cost.rh = fparm(:, 2);
                cost.kk = fparm(:, 3);
                cost.mm = fparm(:, 4);
            end
            dm.userdata.legacy_opf_user_mods = struct( ...
                'z', z, ...
                'lin', lin, ...
                'nlc', {nlc}, ...
                'cost', cost );
        end

        function uc = legacy_user_nln_constraints(obj, dm, mpopt)
            % Handle pre-processing of inputs related to legacy user-defined
            % non-linear constraints, specifically optional ``mpc`` fields
            % ``user_constraints.nle`` and ``user_constraints.nli``.
            %
            % Called by legacy_user_mod_inputs() method.

            mpc = dm.source;

            %% check for user-defined nonlinear constraints
            nnle = 0;   %% number of nonlinear user-defined equality cons
            nnli = 0;   %% number of nonlinear user-defined inequality cons
            if isfield(mpc, 'user_constraints')
                if isfield(mpc.user_constraints, 'nle')
                    for k = 1:length(mpc.user_constraints.nle)
                        nnle = nnle + mpc.user_constraints.nle{k}{2};
                    end
                end
                if isfield(mpc.user_constraints, 'nli')
                    for k = 1:length(mpc.user_constraints.nli)
                        nnli = nnli + mpc.user_constraints.nli{k}{2};
                    end
                end
            end

            %% initialize cell array for add_nln_constraint() args
            uc = cell(nnle+nnli, 1);

            %% user-defined nonlinear equalities
            if nnle
                for k = 1:length(mpc.user_constraints.nle)
                    nlc = mpc.user_constraints.nle{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
                    uc{k} = {nlc{1:2}, 1, fcn, hess, nlc{5}};
                end
            end

            %% user-defined nonlinear inequalities
            if nnli
                for k = 1:length(mpc.user_constraints.nli)
                    nlc = mpc.user_constraints.nli{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})'])
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})'])
                    uc{nnle+k} = {nlc{1:2}, 0, fcn, hess, nlc{5}};
                end
            end
        end
    end     %% methods
end         %% classdef
