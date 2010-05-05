function [results, success, raw] = opf_execute(om, mpopt)
%OPF_EXECUTE  Executes the OPF specified by an OPF model object.
%   [RESULTS, SUCCESS, RAW] = OPF_EXECUTE(OM, MPOPT)
%
%   RESULTS are returned with internal indexing, all equipment
%   in-service, etc.
%
%   See also OPF, OPF_SETUP.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2009-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  setup  -----
%% options
dc  = mpopt(10);        %% PF_DC        : 1 = DC OPF, 0 = AC OPF 
alg = mpopt(11);        %% OPF_ALG
verbose = mpopt(31);    %% VERBOSE

%% build user-defined costs
om = build_cost_params(om);

%% get indexing
[vv, ll, nn] = get_idx(om);

if verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end

%%-----  run DC OPF solver  -----
if dc
  if verbose > 0
    fprintf(' -- DC Optimal Power Flow\n');
  end
  [results, success, raw] = dcopf_solver(om, mpopt);
else
  %%-----  run AC OPF solver  -----
  if verbose > 0
    fprintf(' -- AC Optimal Power Flow\n');
  end

  %% if OPF_ALG not set, choose best available option
  if alg == 0
    if have_fcn('pdipmopf')
      alg = 540;                %% PDIPM
    else
      alg = 560;                %% MIPS
    end
  end

  %% update deprecated algorithm codes to new, generalized formulation equivalents
  if alg == 100 || alg == 200        %% CONSTR
    alg = 300;
  elseif alg == 120 || alg == 220    %% dense LP
    alg = 320;
  elseif alg == 140 || alg == 240    %% sparse (relaxed) LP
    alg = 340;
  elseif alg == 160 || alg == 260    %% sparse (full) LP
    alg = 360;
  end

  %% run specific AC OPF solver
  if alg == 560 || alg == 565                   %% MIPS
    if have_fcn('anon_fcns')
      solver = @mipsopf_solver;
    else
      solver = @mips6opf_solver;
    end
    [results, success, raw] = feval(solver, om, mpopt);
  elseif alg == 540 || alg == 545 || alg == 550 %% PDIPM_OPF, SCPDIPM_OPF, or TRALM_OPF
    if alg == 540                               %% PDIPM_OPF
      if ~have_fcn('pdipmopf')
        error('opf_execute: OPF_ALG %d requires PDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    elseif alg == 545                           %% SCPDIPM_OPF
      if ~have_fcn('scpdipmopf')
        error('opf_execute: OPF_ALG %d requires SCPDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    elseif alg == 550                           %% TRALM_OPF
      if ~have_fcn('tralmopf')
        error('opf_execute: OPF_ALG %d requires TRALM (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    end
    [results, success, raw] = tspopf_solver(om, mpopt);
  elseif alg == 500                             %% MINOPF
    if ~have_fcn('minopf')
      error('opf_execute: OPF_ALG %d requires MINOPF (see http://www.pserc.cornell.edu/minopf/)', alg);
    end
    [results, success, raw] = mopf_solver(om, mpopt);
  elseif alg == 520                             %% FMINCON
    if ~have_fcn('fmincon')
      error('opf_execute: OPF_ALG %d requires FMINCON (Optimization Toolbox 2.x or later)', alg);
    end
    if have_fcn('anon_fcns')
      solver = @fmincopf_solver;
    else
      solver = @fmincopf6_solver;
    end
    [results, success, raw] = feval(solver, om, mpopt);
  elseif alg == 300                             %% CONSTR
    if ~have_fcn('constr')
      error('opf_execute: OPF_ALG %d requires CONSTR (Optimization Toolbox 1.x)', alg);
    end
    [results, success, raw] = copf_solver(om, mpopt);
  elseif alg == 320 || alg == 340 || alg == 360 %% LP
    [results, success, raw] = lpopf_solver(om, mpopt);
  else
    error('opf_execute: OPF_ALG %d is not a valid algorithm code', alg);
  end
end
if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
    raw.output.alg = alg;
end

if success
  if ~dc 
    %% copy bus voltages back to gen matrix
    results.gen(:, VG) = results.bus(results.gen(:, GEN_BUS), VM);

    %% gen PQ capability curve multipliers
    if ll.N.PQh > 0 || ll.N.PQl > 0
      mu_PQh = results.mu.lin.l(ll.i1.PQh:ll.iN.PQh) - results.mu.lin.u(ll.i1.PQh:ll.iN.PQh);
      mu_PQl = results.mu.lin.l(ll.i1.PQl:ll.iN.PQl) - results.mu.lin.u(ll.i1.PQl:ll.iN.PQl);
      Apqdata = userdata(om, 'Apqdata');
      results.gen = update_mupq(results.baseMVA, results.gen, mu_PQh, mu_PQl, Apqdata);
    end

    %% compute g, dg, f, df, d2f if requested by RETURN_RAW_DER = 1
    if mpopt(52)
      %% move from results to raw if using v4.0 of MINOPF or TSPOPF
      if isfield(results, 'dg')
        raw.dg = results.dg;
        raw.g = results.g;
      end
      %% compute g, dg, unless already done by post-v4.0 MINOPF or TSPOPF
      if ~isfield(raw, 'dg')
        mpc = get_mpc(om);
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
        [g, geq, dg, dgeq] = opf_consfcn(results.x, om, Ybus, Yf, Yt, mpopt);
        raw.g = [ geq; g];
        raw.dg = [ dgeq'; dg'];   %% true Jacobian organization
      end
      %% compute df, d2f
      [f, df, d2f] = opf_costfcn(results.x, om);
      raw.df = df;
      raw.d2f = d2f;
    end
  end

  %% delete g and dg fieldsfrom results if using v4.0 of MINOPF or TSPOPF
  if isfield(results, 'dg')
    rmfield(results, 'dg');
    rmfield(results, 'g');
  end

  %% angle limit constraint multipliers
  if ll.N.ang > 0
    iang = userdata(om, 'iang');
    results.branch(iang, MU_ANGMIN) = results.mu.lin.l(ll.i1.ang:ll.iN.ang) * pi/180;
    results.branch(iang, MU_ANGMAX) = results.mu.lin.u(ll.i1.ang:ll.iN.ang) * pi/180;
  end
end

%% assign values and limit shadow prices for variables
om_var_order = get(om, 'var', 'order');
for k = 1:length(om_var_order)
  name = om_var_order{k};
  if getN(om, 'var', name)
    idx = vv.i1.(name):vv.iN.(name);
    results.var.val.(name) = results.x(idx);
    results.var.mu.l.(name) = results.mu.var.l(idx);
    results.var.mu.u.(name) = results.mu.var.u(idx);
  end
end

%% assign shadow prices for linear constraints
om_lin_order = get(om, 'lin', 'order');
for k = 1:length(om_lin_order)
  name = om_lin_order{k};
  if getN(om, 'lin', name)
    idx = ll.i1.(name):ll.iN.(name);
    results.lin.mu.l.(name) = results.mu.lin.l(idx);
    results.lin.mu.u.(name) = results.mu.lin.u(idx);
  end
end

%% assign shadow prices for non-linear constraints
if ~dc
  om_nln_order = get(om, 'nln', 'order');
  for k = 1:length(om_nln_order)
    name = om_nln_order{k};
    if getN(om, 'nln', name)
      idx = nn.i1.(name):nn.iN.(name);
      results.nln.mu.l.(name) = results.mu.nln.l(idx);
      results.nln.mu.u.(name) = results.mu.nln.u(idx);
    end
  end
end

%% assign values for components of user cost
om_cost_order = get(om, 'cost', 'order');
for k = 1:length(om_cost_order)
  name = om_cost_order{k};
  if getN(om, 'cost', name)
    results.cost.(name) = compute_cost(om, results.x, name);
  end
end

%% if single-block PWL costs were converted to POLY, insert dummy y into x
%% Note: The "y" portion of x will be nonsense, but everything should at
%%       least be in the expected locations.
pwl1 = userdata(om, 'pwl1');
if ~isempty(pwl1) && alg ~= 545 && alg ~= 550
  %% get indexing
  vv = get_idx(om);
  if dc
    nx = vv.N.Pg;
  else
    nx = vv.N.Qg;
  end
  y = zeros(length(pwl1), 1);
  raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
  results.x = [ results.x(1:nx); y; results.x(nx+1:end)];
end

