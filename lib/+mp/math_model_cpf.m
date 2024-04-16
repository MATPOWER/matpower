classdef (Abstract) math_model_cpf < mp.math_model_pf
% mp.math_model_cpf - Abstract base class for continuation power flow (CPF) **math model** objects.

%   MATPOWER
%   Copyright (c) 2021-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function tag = task_tag(obj)
            %

            tag = 'cpf';
        end

        function name = task_name(obj)
            %

            name = 'Continuation Power Flow';
        end

        function obj = add_aux_data(obj, nm, dm, mpopt)
            %

            %% create aux_data struct
            obj.aux_data = obj.build_aux_data_cpf(nm, dm, mpopt);
        end

        %% can't simply override build_aux_data(), since we
        %% need to be able to call the PF version too
        function ad = build_aux_data_cpf(obj, nm, dm, mpopt)
            %

            dmt = dm.userdata.target;
            nmt = nm.userdata.target;

            ad  = obj.build_aux_data(nm, dm, mpopt);
            adt = obj.build_aux_data(nmt, dmt, mpopt);

            ad = obj.check_xfer(nm, nmt, ad, adt);
            ad.nmt = nmt;
            ad.dmb = dm.copy(); %% save copy of unmodified base case data model
            ad.dmt = dmt;
            ad.adt = adt;
        end

        function add_vars(obj, nm, dm, mpopt)
            %

            add_vars@mp.math_model_pf(obj, nm, dm, mpopt);
            %% required for using nmt to compute nm state from mm state
            obj.aux_data.adt.var_map = obj.aux_data.var_map;
            obj.add_var('lambda', 1, 0);
        end

        %% can't simply override convert_x_m2n(), since we
        %% need to be able to call the PF version too
        function [vx_, z_, x_] = convert_x_m2n_cpf(obj, mmx, nm, only_v)
            %

            ad = obj.aux_data;
            nmt = ad.nmt;
            lam = mmx(end);     %% continuation parameter lambda

            %% update voltages and get base z_
            [vx_,  zb_] = obj.convert_x_m2n(mmx(1:end-1), nm, 1);
            obj.aux_data = ad.adt;  %% swap to target aux_data
            [vxt_, zt_] = obj.convert_x_m2n(mmx(1:end-1), nmt, 1);
            obj.aux_data = ad;      %% restore base aux_data
            assert(norm(vx_-vxt_, Inf) < eps);

            %% compute z_ as function of continuation parameter lambda
            z_ = (1-lam) * zb_+ lam * zt_;

            %% update dependent portions of z, if requested
            if nargin < 4 || ~only_v
                rpv = [ad.ref; ad.pv];      %% slack and PV nodes
                idx = find(any(nm.C(rpv, :), 1));   %% ports connected to slack/PV nodes
                Sinjb =  nm.port_inj_power([vx_; zb_], 1, idx);
                Sinjt = nmt.port_inj_power([vx_; zt_], 1, idx);
                Sinj = (1-lam) * Sinjb + lam * Sinjt;
                z_ = obj.update_z(nm, vx_, z_, ad, Sinj, idx);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        %% can't simply override node_balance_equations(), since we
        %% need to be able to call the PF version too
        function [f, J] = node_balance_equations_cpf(obj, x, nm)
            %

            ad = obj.aux_data;
            nmt = ad.nmt;
            lam = x(end);   %% continuation parameter lambda

            if nargout > 1
                [fb, Jb] = obj.node_balance_equations(x(1:end-1), nm);
                obj.aux_data = ad.adt;  %% swap to target aux_data
                [ft, Jt] = obj.node_balance_equations(x(1:end-1), nmt);
                obj.aux_data = ad;      %% restore base aux_data
                J = [(1-lam) * Jb + lam * Jt    ft - fb];
            else
                fb = obj.node_balance_equations(x(1:end-1), nm);
                obj.aux_data = ad.adt;  %% swap to target aux_data
                ft = obj.node_balance_equations(x(1:end-1), nmt);
                obj.aux_data = ad;      %% restore base aux_data
            end
            f = (1-lam) * fb + lam * ft;
        end

        function nm = network_model_x_soln(obj, nm)
            %

            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                obj.convert_x_m2n_cpf(obj.soln.x, nm);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %

            ad = obj.aux_data;
            opt = mpopt2pneopt(mpopt);
            opt.output_fcn = @(varargin)obj.pne_output_fcn(nm, ad, varargin{:});
            opt.plot.idx_default = @()obj.plot_idx_default(nm, dm, ad);
            opt.plot.yfcn = @(v_,idx)obj.plot_yfcn(nm, dm, ad, v_, idx);
            opt = obj.add_callbacks(opt, nm, dm, mpopt);
        end

        function ad = check_xfer(obj, nm, nmt, ad, adt)
            %

            %% ensure base and target parameters are identical,
            %% except for fixed power injections & states
            [L,  N,  i,  s ] = nm.get_params([], {'L', 'N', 'i', 's'});
            [Lt, Nt, it, st] = nmt.get_params([], {'L', 'N', 'i', 's'});
            tol = 1e-12;

            %% nodal transfer from direct power injection states
            NN = nm.C * N;
            NNt = nm.C * Nt;
            LL = nm.C * L;
            LLt = nm.C * Lt;
            zzr = ad.zr - adt.zr;
            zzi = ad.zi - adt.zi;
            kref = find(any(LL(ad.ref, :), 1) | any(NN(ad.ref, :), 1));
            kpv  = find(any(LL(ad.pv,  :), 1) | any(NN(ad.pv,  :), 1));
            zzr(kref) = 0;   %% zero active transfer at slack node
            zzi(kpv)  = 0;   %% zero reactive transfer at PV nodes
            zz = zzr + 1j * zzi;

            %% nodal transfer from constant power injections
            ss = nm.C * (s - st);

            %% create transfer vector from diff in direct power injections
            %% between base and target, from constant power elements and
            %% direct power injection states, used only to auto select
            %% largest transfer bus for default voltage nose-curve plot
            ad.xfer = ss + NN * zz;

            %% Power flow equations must be linear in continuation parameter
            %% lambda. To do that we must ensure that ...
            %% 1. Specified voltages do not vary with lambda.
            [va,  vm ] = nm.aux_data_va_vm(ad);
            [vat, vmt] = nm.aux_data_va_vm(adt);
            rpv = [ad.ref; ad.pv];
            if norm(va(ad.ref)-vat(ad.ref), Inf) > tol
                error('mp.math_model_cpf.check_xfer: base and target cases must have identical voltage angles at reference nodes.')
            end
            if norm(vm(rpv)-vmt(rpv), Inf) > tol
                error('mp.math_model_cpf.check_xfer: base and target cases must have identical voltage magnitudes at reference and PV nodes.')
            end
            %% 2. Elements of z that vary with lambda must have only constant
            %%    coefficients, i.e. corresponding columns of L and N must be
            %%    identical in base and target models.
            k = find(zz);
            if norm(LL(:,k) - LLt(:,k), Inf) > tol
                error('mp.math_model_cpf.check_xfer: base and target cases must have identical coefficients for any current injection state variables that vary from base to target.')
            end
            if norm(NN(:,k) - NNt(:,k), Inf) > tol
                error('mp.math_model_cpf.check_xfer: base and target cases must have identical coefficients for any power injection state variables that vary from base to target.')
            end
        end

        function [names, vals] = pne_output_fcn(obj, nm, ad, x, x_hat)
            %

            %% [names, vals] = obj.pne_output_fcn(nm, ad, x, x_hat)
            %% names = obj.pne_output_fcn(nm, ad)
            names = {'V_hat', 'V'};
            if nargin > 3
                [V_hat, ~] = obj.convert_x_m2n_cpf(x_hat, nm, 1);
                [V,     ~] = obj.convert_x_m2n_cpf(x,     nm, 1);
                vals = {V_hat, V};
            end
        end

        function y = plot_yfcn(obj, nm, dm, ad, v_, bus_num)
            %

            %% find node idx from external bus number
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping
            nidx = nm.get_node_idx('bus');

            k = find( bus_num < 0 | bus_num > length(b2i) );
            if ~isempty(k)
                error('mp.math_model_cpf.plot_yfcn: %d is not a valid bus number for CPF voltage plot', bus_idx(k));
            end

            idx = nidx(b2i(bus_num));
            y = abs(v_(idx, :));
        end

        function idx = plot_idx_default(obj, nm, dm, ad)
            %

            %% plot voltage of PQ node with max transfer as default
            nidx = nm.get_node_idx('bus');     %% node indices of buses
            [~, i] = max(abs(ad.xfer(ad.pq)) .* ismember(ad.pq, nidx));
            bi = ad.pq(i);                      %% index of bus w/max transfer
            idx = dm.elements.bus.ID(bi);       %% bus num of same bus
        end

        function opt = add_callbacks(obj, opt, nm, dm, mpopt)
            %

            qlim = mpopt.cpf.enforce_q_lims;    %% enforce reactive limits
            plim = mpopt.cpf.enforce_p_lims;    %% enforce active limits
            vlim = mpopt.cpf.enforce_v_lims;    %% enforce voltage magnitude limits
            flim = mpopt.cpf.enforce_flow_lims; %% enforce branch flow limits

            %% initialize event and callback options
            if ~isfield(opt, 'events') || isempty(opt.events)
                opt.events = {};
            end
            if ~isfield(opt, 'callbacks') || isempty(opt.callbacks)
                opt.callbacks = {};
            end

            if flim
                opt.events{end+1} = { ...
                    'FLIM', ...
                    @(cx, opt)event_flim(obj, cx, opt, nm, dm, mpopt), ...
                    mpopt.cpf.flow_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)callback_flim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt), ...
                    53 };
            end

            if vlim
                opt.events{end+1} = { ...
                    'VLIM', ...
                    @(cx, opt)event_vlim(obj, cx, opt, nm, dm, mpopt), ...
                    mpopt.cpf.v_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)callback_vlim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt), ...
                    52 };
            end

            if qlim
                opt.events{end+1} = { ...
                    'QLIM', ...
                    @(cx, opt)event_qlim(obj, cx, opt, nm, dm, mpopt), ...
                    mpopt.cpf.q_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)callback_qlim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt), ...
                    41 };
            end

            if plim
                opt.events{end+1} = { ...
                    'PLIM', ...
                    @(cx, opt)event_plim(obj, cx, opt, nm, dm, mpopt), ...
                    mpopt.cpf.p_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)callback_plim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt), ...
                    40 };
            end
        end

        function efv = event_flim(obj, cx, opt, nm, dm, mpopt)
            %

            %% get branch flow constraints
            branch_nme = nm.elements.branch;
            branch_dme = branch_nme.data_model_element(dm);
            rate_a = branch_dme.rate_a * dm.base_mva;
            ibr = find(rate_a ~= 0 & rate_a < 1e10);
            nl2 = length(ibr);          %% number of constrained branches

            if nl2
                nl = branch_nme.nk;     %% port indexes

                %% convert cx.x back to x_
                x_ = obj.convert_x_m2n_cpf(cx.x, nm);

                %% branch flows
                S_fr = branch_nme.port_inj_power(x_, 1, ibr)    * dm.base_mva;
                S_to = branch_nme.port_inj_power(x_, 1, nl+ibr) * dm.base_mva;
                S_fr = sqrt(S_fr .* conj(S_fr));
                S_to = sqrt(S_to .* conj(S_to));

                %% branch flow lim event function
                efv = max(S_fr, S_to) - rate_a(ibr);
            else
                efv = NaN;
            end
        end

        function efv = event_qlim(obj, cx, opt, nm, dm, mpopt)
            %

            ad = obj.aux_data;

            %% convert cx.x back to v_, z_
            [v_, z_] = obj.convert_x_m2n_cpf(cx.x, nm);

            %% coefficient matrix for power injection states
            rpv = [ad.ref; ad.pv];      %% slack and PV nodes
            CCrpv = nm.C(rpv, :) * nm.get_params([], 'N') * nm.D';
            jrpv = find(any(CCrpv, 1)); %% indices for states at slack/PV nodes

            %% limit violations at slack/PV nodes
            [~, zi_min, zi_max] = nm.params_var('zi'); %% bounds on zi
            v_Qmax = NaN(nm.nz, 1);
            v_Qmin = v_Qmax;
            v_Qmax(jrpv) = imag(z_(jrpv)) - zi_max(jrpv);
            v_Qmin(jrpv) = zi_min(jrpv) - imag(z_(jrpv));

            %% assemble event function value
            efv = [v_Qmax; v_Qmin] * dm.base_mva;
        end

        function efv = event_plim(obj, cx, opt, nm, dm, mpopt)
            %

            %% convert cx.x back to v_, z_
            [v_, z_] = obj.convert_x_m2n_cpf(cx.x, nm);

            %% limit violations
            [~, ~, zr_max] = nm.params_var('zr'); %% bounds on zr
            v_Pmax = real(z_) - zr_max;

            %% ignore those that are already at their max limit
            if isfield(cx.cbs, 'plim') && ~isempty(cx.cbs.plim.idx)
                v_Pmax(cx.cbs.plim.idx) = NaN;
            end

            %% assemble event function value
            efv = v_Pmax * dm.base_mva;
        end

        function [nx, cx, s] = callback_flim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt)
            %

            %% initialize
            if k == 0   %% check for base case flow violations
                %% get branch flow constraints
                branch_nme = nm.elements.branch;
                branch_dme = branch_nme.data_model_element(dm);
                rate_a = branch_dme.rate_a * dm.base_mva;
                ibr = find(rate_a ~= 0 & rate_a < 1e10);
                nl2 = length(ibr);          %% number of constrained branches

                if nl2
                    nl = branch_nme.nk;     %% port indexes

                    %% convert cx.x back to x_
                    x_ = obj.convert_x_m2n_cpf(cx.x, nm);

                    %% branch flows
                    S_fr = branch_nme.port_inj_power(x_, 1, ibr)    * dm.base_mva;
                    S_to = branch_nme.port_inj_power(x_, 1, nl+ibr) * dm.base_mva;
                    S_fr = sqrt(S_fr .* conj(S_fr));
                    S_to = sqrt(S_to .* conj(S_to));

                    %% violated branch flows
                    if any(max(S_fr, S_to) > rate_a(ibr))
                        %% find the lines and which lim(s)
                        iL = find(max(S_fr, S_to) > rate_a(ibr));
                        msg = '';
                        for j = 1:length(iL)
                            L = ibr(iL(j));
                            fidx = find(branch_nme.C(:, L));
                            tidx = find(branch_nme.C(:, nl+L));
                            flabel = nm.set_type_label('node', fidx, dm);
                            tlabel = nm.set_type_label('node', tidx, dm);

                            msg = sprintf('%sbranch flow limit violated in base case: %s -- %s exceeds limit of %g MVA\n',...
                                msg, flabel, tlabel, rate_a(L));
                        end

                        %% prepare to terminate
                        s.done = 1;
                        s.done_msg = msg;
                    end
                end
            end

            %% skip if finalize or done
            if k < 0 || s.done
                return;
            end

            %% handle event
            ev = pne_detected_event(s.events, 'FLIM', 1);   %% zero only
            if ~isempty(ev)
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end

                %% get branch flow constraints
                branch_nme = nm.elements.branch;
                branch_dme = branch_nme.data_model_element(dm);
                rate_a = branch_dme.rate_a * dm.base_mva;
                ibr = find(rate_a ~= 0 & rate_a < 1e10);
                nl2 = length(ibr);      %% number of constrained branches
                nl = branch_nme.nk;     %% port indexes

                %% find branch(es) with violated lim(s)
                iL = ev.idx;            %% event function index
                for j = 1:length(iL)
                    L = ibr(iL(j)); %% index of critical branch event of interest
                    fidx = find(branch_nme.C(:, L));
                    tidx = find(branch_nme.C(:, nl+L));
                    flabel = nm.set_type_label('node', fidx, dm);
                    tlabel = nm.set_type_label('node', tidx, dm);

                    msg = sprintf('%sbranch flow limit reached\nbranch: %s -- %s at limit of %g MVA @ lambda = %.4g, in %d continuation steps',...
                        msg, flabel, tlabel, rate_a(L), nx.x(end), k);
                end

                %% prepare to terminate
                s.done = 1;
                s.done_msg = msg;
            end
        end

        function [nx, cx, s] = callback_qlim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt)
            %

            %% skip if initialize, finalize or done
            if k <= 0 || s.done
                return;
            end

            %% handle event
            [ev, i] = pne_detected_event(s.events, 'QLIM', 1);  %% zero only
            if ~isempty(ev)
                ad = obj.aux_data;
                if ad.nref ~= 1
                    error('mp.math_model_cpf.callback_qlim: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus');
                end

                efidx = ev.idx;             %% event function index
                [~, zi_min, zi_max] = nm.params_var('zi');  %% bounds on zi
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end
                for j = 1:length(efidx)
                    %% find index of z var
                    idx = efidx(j);         %% index of z var
                    if idx <= nm.nz
                        maxlim = 1;         %% Qmax violation
                        lim = zi_max(idx) * dm.base_mva;
                        lim_type = 'Qmax';
                    else
                        idx = idx - nm.nz;  %% correct index of z var
                        maxlim = 0;         %% Qmin violation
                        lim = zi_min(idx) * dm.base_mva;
                        lim_type = 'Qmin';
                    end

                    %% get label for z var
                    zlabel = nm.set_type_label('state', idx, dm);

                    %% get label for corresponding node
                    CC = nm.C * nm.get_params([], 'N') * nm.D';
                    nidx = find(CC(:, idx));
                    nlabel = nm.set_type_label('node', nidx, dm);

                    msg = sprintf('%s%s @ %s reached %g MVAr %s lim @ lambda = %.4g : %s converted to PQ', ...
                            msg, zlabel, nlabel, lim, lim_type, ...
                            nx.x(end), nlabel);

                    %% set Q to exact limit
                    [v_, z_] = obj.convert_x_m2n_cpf(nx.x, nm);
                    z_(idx) = real(z_(idx)) + 1j * lim / dm.base_mva;

                    %% change node type to PQ
                    nm.set_node_type_pq(dm, nidx);

                    %% check for existence of remaining slack/PV bus
                    try
                        %% potentially pick new reference bus
                        [ref, pv, pq] = nm.node_types(nm, dm);
                    catch
                        s.done = 1;
                        s.done_msg = 'No REF or PV nodes remaining.';

                        %% undo change of last REF to PQ
                        nm.set_node_type_ref(dm, ad.ref);

                        break;
                    end

                    if ~s.done
                        %% get target case data, network models
                        dmt = ad.dmt;
                        nmt = ad.nmt;

                        %% change node type in target case
                        nmt.set_node_type_pq(dmt, nidx);

                        %% zero out Q transfer for bus
                        ss = nm.set_type_idx_map('zi', idx);
                        nm.zi.data.v0.(ss.name)(ss.i) = imag(z_(idx));
                        nmt.zi.data.v0.(ss.name)(ss.i) = imag(z_(idx));

                        %% if slack changed ...
                        if ref ~= ad.ref
                            %% find zr corresponding to all zr at ref node
                            zref = find(CC(ad.ref, :));
                            ss = nm.set_type_idx_map('zr', zref);

                            %% zero out P transfer at old ref
                            for kk = 1:length(ss)
                                nm.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(zref(kk)));
                                nmt.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(zref(kk)));
                            end

                            %% update voltage angle at new ref node
                            ss = nm.set_type_idx_map('va', ref);
                            nm.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                            nmt.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                        end
                    end
                end
                if ~s.done
                    s.done = 1;
                    dir_from_jac_eigs = isempty(find(nx.z(ad.npv+ad.npq+1:end-1) > 0, 1));
                    s.warmstart = struct('nmt', nmt, 'dmt', dmt, ...
                        'dir_from_jac_eigs', dir_from_jac_eigs);
                end
                s.events(i).msg = msg;
            end
        end

        function [nx, cx, s] = callback_plim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt)
            %

            %% skip if finalize or done
            if k < 0 || s.done
                return;
            elseif k == 0 && ~isfield(cx.cbs, 'plim')
                cx.cbs.plim.idx = [];
            end

            %% handle event
            [ev, i] = pne_detected_event(s.events, 'PLIM', 1);  %% zero only
            if ~isempty(ev)
                ad = obj.aux_data;
                if ad.nref ~= 1
                    error('mp.math_model_cpf.callback_plim: ''cpf.enforce_plims'' option only valid for systems with exactly one REF bus');
                end

                efidx = ev.idx;             %% event function index
                [~, ~, zr_max] = nm.params_var('zr'); %% bounds on zr
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end
                for j = 1:length(efidx)
                    %% find index of z var
                    idx = efidx(j);         %% index of z var
                    lim = zr_max(idx) * dm.base_mva;

                    %% get label for z var
                    zlabel = nm.set_type_label('state', idx, dm);

                    %% get label for corresponding node
                    CC = nm.C * nm.get_params([], 'N') * nm.D';
                    nidx = find(CC(:, idx));
                    nlabel = nm.set_type_label('node', nidx, dm);

                    msg = sprintf('%s%s @ %s reached %g MW pg upper bound @ lambda = %.4g', ...
                            msg, zlabel, nlabel, lim, nx.x(end));

                    %% set P to exact limit
                    [v_, z_] = obj.convert_x_m2n_cpf(nx.x, nm);
                    z_(idx) = lim / dm.base_mva + 1j * imag(z_(idx));

                    dmt = ad.dmt;
                    nmt = ad.nmt;

                    %% find zr corresponding to all zr at this node
                    k = find(CC(nidx, :));
                    ss = nm.set_type_idx_map('zr', k);

                    %% set z to limit at this node
                    for kk = 1:length(ss)
                        nm.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(k(kk)));
                        nmt.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(k(kk)));
                    end

                    %% save index of element at limit
                    nx.cbs.plim.idx = [nx.cbs.plim.idx; idx];

                    %% if it is at the ref node
                    if nidx == ad.ref
                        %% find free injections (not at max lim)
                        free = ones(nm.nz, 1);
                        free(nx.cbs.plim.idx) = 0;

                        %% first node with any free injection
                        ref = find(any(CC * spdiags(free, 0, nm.nz, nm.nz), 2));
                        if isempty(ref)
                            s.done = 1;
                            s.done_msg = 'All generators at Pmax';
                            break;
                        else
                            %% convert this node to PV, new ref bus to REF
                            nm.set_node_type_pv(dm, nidx);
                            nmt.set_node_type_pv(dmt, nidx);
                            nm.set_node_type_ref(dm, ref);
                            nmt.set_node_type_ref(dmt, ref);

                            %% update voltage angle at new ref node
                            ss = nm.set_type_idx_map('va', ref);
                            nm.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                            nmt.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));

                            rlabel = nm.set_type_label('node', ref, dm);
                            msg = sprintf('%s : ref changed from %s to %s', ...
                                msg, nlabel, rlabel);
                        end
                    end
                end
                if ~s.done
                    s.done = 1;
                    s.warmstart = struct('nmt', nmt, 'dmt', dmt);
                end
                s.events(i).msg = msg;
            end
        end
    end     %% methods
end         %% classdef
