function hh = plot_storage(md, idx, varargin)
%PLOT_STORAGE   Plot storage unit results
%
%   PLOT_STORAGE(MD)
%   PLOT_STORAGE(MD, IDX)
%   PLOT_STORAGE(MD, IDX, '<option1_name>', '<option1_value', ...)
%   PLOT_STORAGE(MD, IDX, OPT)
%   PLOT_STORAGE(MD, IDX, OPT, '<option1_name>', '<option1_value', ...)
%   H = PLOT_STORAGE(MD, ...)
%
%   IDX is the storage unit index. If IDX is a vector, it sums them first,
%   if empty, it includes all storage units in MD. Options can be
%   specified as an OPT struct or as individual pairs of 'name' and 'value'
%   arguments. The possible options include the following, where the default
%   is shown in parenthesis:
%       'saveit'        (false) flag to indicate whether to create PDF file
%       'saveall'       (false) flag to indicate whether to create individual
%                       PDF files for each element of IDX, as well as
%                       aggregate, when IDX is a vector (or empty)
%       'savepath'      ('') path to directory to save files in
%       'savename'      ('stored-energy-%s.pdf') name of PDF file
%                       %s is optional placeholder for storage unit index
%       'separation'    (0.8)  %% separation of beginning/end of period (0-1)
%           0   = beginning & end of period t, staircase (dispatches are
%                 vertical lines at t)
%           0.5 = evenly separated (both dispatch and transitions visible)
%           1   = beginning of t aligned with end of t-1, smooth,
%                 (transitions are vertical displacements at t+/-0.5)
%       'sort_tol'      (1e-6) round to nearest sort_tol for sorting
%       'size_factor'   (1) to scale font/marker sizes in case you want to
%                       do sub-plots or something
%       'show_grid'                 (true) vertical lines to divide periods
%       'show_expected_initial'     (true)
%       'show_expected_final'       (true)
%       'show_adjusted_dispatches'  (true)
%       'show_dispatches'           (false)
%
%   Returns handle to current figure window.

%   TO DO: Do initial invisible plot to get v axis parameters.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% initialize some variables
my_xlabel = 'Period';
nt = md.idx.nt;
ns = md.idx.ns;
nj_max = max(md.idx.nj);
baseMVA = md.mpc.baseMVA;

%% input args
if nargin < 2
    idx = [];
end
if isempty(idx)
    idx = (1:ns)';
end
nidx = length(idx);
if nidx > 1 && size(idx, 1) == 1
    idx = idx';     %% convert row vector to column vector
end
g = md.Storage.UnitIdx(idx);
b = md.mpc.gen(g, GEN_BUS);

%% default options
opt = struct( ...
    'saveit', false, ...
    'saveall', false, ...
    'savepath', '', ...
    'savename', 'stored-energy-%s.pdf', ...
    'separation', 0.8, ...
    'sort_tol', 1e-6, ...
    'size_factor', 1, ...
    'show_grid', true, ...
    'show_adjusted_dispatches', true, ...
    'show_expected_initial', true, ...
    'show_expected_final', true, ...
    'show_dispatches', false );

%% process options
if mod(length(varargin), 2) %% odd number of options, first must be OPT struct
    if ~isstruct(varargin{1})
        error('plot_storage: Single OPT argument must be a struct');
    end
    myopt = varargin{1};
    k = 2;
else                        %% even number of options
    myopt = struct;
    k = 1;
end
while k < length(varargin)
    opt_name = varargin{k};
    opt_val  = varargin{k+1};
    if ~isfield(opt, opt_name)
        error('plot_name: ''%s'' is not a valid option name', opt_name);
    end
    myopt.(opt_name) = opt_val;
    k = k + 2;
end
fields = fieldnames(myopt);
for f = 1:length(fields)
    opt.(fields{f}) = myopt.(fields{f});
end

%% call recursively for individual plots if indicated
if opt.saveall && nidx > 1
    if isempty(strfind(opt.savename, '%s'))
        error('plot_storage: ''savename'' must include a ''%%s'' placeholder when ''saveall'' option is true.');
    end
    for i = 1:nidx
        plot_storage(md, idx(i), opt, 'saveit', true);
    end
end

%% extract and expand storage parameters
MinStorageLevel = md.Storage.MinStorageLevel;
MaxStorageLevel = md.Storage.MaxStorageLevel;
if size(MinStorageLevel, 1) == 1 && ns > 1  %% expand rows
  MinStorageLevel = ones(ns, 1) * MinStorageLevel;
end
if size(MinStorageLevel, 2) == 1 && nt > 1  %% expand cols
  MinStorageLevel = MinStorageLevel * ones(1, nt);
end
if size(MaxStorageLevel, 1) == 1 && ns > 1  %% expand rows
  MaxStorageLevel = ones(ns, 1) * MaxStorageLevel;
end
if size(MaxStorageLevel, 2) == 1 && nt > 1  %% expand cols
  MaxStorageLevel = MaxStorageLevel * ones(1, nt);
end
if isempty(md.Storage.InEff)
  InEff = 1;                        %% no efficiency loss by default
else
  InEff = md.Storage.InEff;
end
if size(InEff, 1) == 1 && ns > 1    %% expand rows
  InEff = ones(ns, 1) * InEff;
end
if size(InEff, 2) == 1 && nt > 1    %% expand cols
  InEff = InEff * ones(1, nt);
end
if isempty(md.Storage.OutEff)
  OutEff = 1;                       %% no efficiency loss by default
else
  OutEff = md.Storage.OutEff;
end
if size(OutEff, 1) == 1 && ns > 1   %% expand rows
  OutEff = ones(ns, 1) * OutEff;
end
if size(OutEff, 2) == 1 && nt > 1   %% expand cols
  OutEff = OutEff * ones(1, nt);
end
if isempty(md.Storage.LossFactor)
  LossFactor = 0;                       %% no losses by default
else
  LossFactor = md.Storage.LossFactor;
end
if size(LossFactor, 1) == 1 && ns > 1   %% expand rows
  LossFactor = ones(ns, 1) * LossFactor;
end
if size(LossFactor, 2) == 1 && nt > 1   %% expand cols
  LossFactor = LossFactor * ones(1, nt);
end
if isempty(md.Storage.rho)
  rho = 1;                      %% use worst case by default (for backward compatibility)
else
  rho = md.Storage.rho;
end
if size(rho, 1) == 1 && ns > 1  %% expand rows
  rho = ones(ns, 1) * rho;
end
if size(rho, 2) == 1 && nt > 1  %% expand cols
  rho = rho * ones(1, nt);
end

%% initialize data structures to be plotted
offset = opt.separation/2;          %% offset for plotting on horiz axis
p = (1:nt)';
pp = 2*p;
ppp(pp-1) = p - offset;
ppp(pp  ) = p + offset;
Sp = md.results.Sp(idx, :);                     %% s+
Sm = md.results.Sm(idx, :);                     %% s-
eS = md.Storage.ExpectedStorageState(idx, :);   %% expected s
Smin = MinStorageLevel(idx, :);     %% physical lower bound
Smax = MaxStorageLevel(idx, :);     %% physical upper bound
SSp   = NaN(2*nt, 1);               %% expanded for 2 pts per period
SSm   = NaN(2*nt, 1);               %%    "
eeS   = NaN(2*nt, 1);               %%    "
SSmin(pp)   = sum(Smin, 1);
SSmax(pp)   = sum(Smax, 1);
SSmin(pp-1) = SSmin(pp);
SSmax(pp-1) = SSmax(pp);

eSi = NaN(nj_max, nt, nidx);
eSf = NaN(nj_max, nt, nidx);
Si_min = NaN(nj_max, nt, nidx);
Si_max = NaN(nj_max, nt, nidx);
dSi = zeros(nj_max, nt, nidx);  %% deltaS, total change in stored energy from injections
Pg = zeros(nj_max, nt, nidx);   %% total dispatch of specified storage units
jmin = zeros(1, nt);
jmax = zeros(1, nt);
LossCoeff = md.Delta_T * LossFactor/2;
beta1 = (1-LossCoeff) ./ (1+LossCoeff);
beta2 = 1 ./ (1+LossCoeff);
for t = 1:nt
    if t == 1
        if md.Storage.ForceCyclicStorage
            SSm(t) = sum(md.Storage.InitialStorage(idx));
            SSp(t) = sum(md.Storage.InitialStorage(idx));
        else
            SSm(t) = sum(md.Storage.InitialStorageLowerBound(idx));
            SSp(t) = sum(md.Storage.InitialStorageUpperBound(idx));
        end
        eeS(t) = sum(md.Storage.InitialStorage(idx));
    else
        SSm(2*t-1) = sum(Sm(:, t-1), 1);
        SSp(2*t-1) = sum(Sp(:, t-1), 1);
        eeS(2*t-1) = sum(eS(:, t-1), 1);
    end
    %% end of period
    SSm(2*t)   = sum(Sm(:, t), 1);
    SSp(2*t)   = sum(Sp(:, t), 1);
    eeS(2*t)   = sum(eS(:, t), 1);

    vv = get_idx(md.om);
    for j = 1:md.idx.nj(t)
        dS = -md.Delta_T * ...
            (InEff(:,t)  .* md.QP.x(vv.i1.Psc(t,j,1)-1+(1:ns)) + ...
             1./OutEff(:,t) .* md.QP.x(vv.i1.Psd(t,j,1)-1+(1:ns)) );
        dSi(j,t,:) = dS(idx) * baseMVA;
        Pg(j,t,:)  = md.flow(t,j,1).mpc.gen(g, PG);
        Lij = md.tstep(t).Li( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :);
        Lfj = md.tstep(t).Lf( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :);
        Mj  = md.tstep(t).Mg( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :) + ...
              md.tstep(t).Mh( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :);
        Nj  = md.tstep(t).Ng( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :) + ...
              md.tstep(t).Nh( j:md.idx.nj(t):(ns-1)*md.idx.nj(t)+j, :);
        temp_eSi = Lij * md.Storage.InitialStorage/baseMVA + Mj * md.QP.x;
        temp_eSf = Lfj * md.Storage.InitialStorage/baseMVA + Nj * md.QP.x;
        eSi(j,t,:) = baseMVA * temp_eSi(idx);
        eSf(j,t,:) = baseMVA * temp_eSf(idx);
        if t == 1
            if md.Storage.ForceCyclicStorage
                Si_min(j,t,:) = md.Storage.InitialStorage(idx);
                Si_max(j,t,:) = md.Storage.InitialStorage(idx);
            else
                Si_min(j,t,:) = md.Storage.InitialStorageLowerBound(idx);
                Si_max(j,t,:) = md.Storage.InitialStorageUpperBound(idx);
            end
        else
            tmp_eSi = reshape(eSi(j,t,:), nidx, 1);
            Si_min(j,t,:) = rho(idx,t) .* Sm(:,t-1) + (1-rho(idx,t)) .* tmp_eSi;
            Si_max(j,t,:) = rho(idx,t) .* Sp(:,t-1) + (1-rho(idx,t)) .* tmp_eSi;
        end
    end
%     if nidx == -1
%         rows = opt.sort_tol * round(opt.sort_tol^-1 * [ beta1(idx,t)*Si_min(:,t,1) + beta2(idx,t)*dSi(:,t,1) Si_min(:,t,1) ]);
%         [junk, iii] = sortrows(rows, [1 2]);
%         jmin(t) = iii(1);
%         rows = opt.sort_tol * round(opt.sort_tol^-1 * [ beta1(idx,t)*Si_max(:,t,1) + beta2(idx,t)*dSi(:,t,1) Si_max(:,t,1) ]);
%         [junk, iii] = sortrows(rows, [-1 -2]);
%         jmax(t) = iii(1);
%     else
        tmpS = reshape(Si_min(:,t,:), nj_max, nidx);
        rows = opt.sort_tol * round(opt.sort_tol^-1 * [tmpS * beta1(idx,t) + reshape(dSi(:,t,:), nj_max, nidx) * beta2(idx,t) tmpS]);
        [junk, i] = sortrows(rows, [1 2]);
        jmin(t) = i(1);
        tmpS = reshape(Si_max(:,t,:), nj_max, nidx);
        rows = opt.sort_tol * round(opt.sort_tol^-1 * [tmpS * beta1(idx,t) + reshape(dSi(:,t,:), nj_max, nidx) * beta2(idx,t) tmpS]);
        [junk, i] = sortrows(rows, [-1 -2]);
        jmax(t) = i(1);
%     end
end
SSi_min = zeros(nt,1);
SSi_max = zeros(nt,1);
for t = 1:nt
    SSi_min(t) = sum(Si_min(jmin(t), t, :), 3);
    SSi_max(t) = sum(Si_max(jmax(t), t, :), 3);
end

%% do plots
%% figure out axis limits
plot(ppp, SSmax+0.001, 'LineStyle', 'none')
hold on
plot(ppp, SSmin-0.001, 'LineStyle', 'none')

%% draw grid
v = axis;
if opt.show_grid
    m = v(3) - 100*(v(4) - v(3));
    M = v(4) + 100*(v(4) - v(3));
    for k = 0:nt
        line([k; k], [m; M], 'LineWidth', 0.25, 'Color', 0.9*[1 1 1]);
    end
end
axis(v);

plot(ppp, SSmin, ':k', 'LineWidth', 1);
plot(ppp, SSmax, ':k', 'LineWidth', 1);
for t = 1:nt
    for j = 1:md.idx.nj(t)
        if opt.show_adjusted_dispatches
            Si = sum(Si_min(j,t,:), 3);
            Sf = reshape(Si_min(j,t,:), 1, nidx) * beta1(idx, t) + ...
                    reshape(dSi(j,t,:), 1, nidx) * beta2(idx,t);
            line(t+offset*[-1;1], [Si; Sf], 'Color', 0.8*[1 1 1]);
            Si = sum(Si_max(j,t,:), 3);
            Sf = reshape(Si_max(j,t,:), 1, nidx) * beta1(idx, t) + ...
                    reshape(dSi(j,t,:), 1, nidx) * beta2(idx,t);
            line(t+offset*[-1;1], [Si; Sf], 'Color', 0.8*[1 1 1]);
        end

        if opt.show_dispatches
            Si = sum(Si_min(j,t,:), 3);
            Sf = Si - sum(Pg(j,t,:), 3);
            line(t+offset*[-1;1], [Si; Sf], 'Color', 0.8*[1 0.7 1]);
            Si = sum(Si_max(j,t,:), 3);
            Sf = Si - sum(Pg(j,t,:), 3);
            line(t+offset*[-1;1], [Si; Sf], 'Color', 0.8*[1 0.7 1]);
        end
    end
end
plot(ppp, SSp, '--', 'LineWidth', 1, 'Color', 0.8*[0 1 0]);
plot(ppp, SSm, '--r', 'LineWidth', 1, 'Color', 0.9*[1,0,0]);
plot(p+offset, SSp(pp), 'v', 'LineWidth', 1, 'Color', 0.8*[0 1 0], 'MarkerSize', 6*opt.size_factor);
plot(p+offset, SSm(pp), '^r', 'LineWidth', 1, 'Color', 0.9*[1,0,0], 'MarkerSize', 6*opt.size_factor);
plot(ppp, eeS, 'Color', [0 0 1], 'LineWidth', 2);
if opt.show_expected_initial
    if any(any(rho(idx,:) > 0))
        plot(p-offset, sum(Si_min, 3), '+', 'LineWidth', 1, 'Color', 0.9*[1,0,0], 'MarkerSize', 4*opt.size_factor);
        plot(p-offset, sum(Si_max, 3), '+', 'LineWidth', 1, 'Color', 0.8*[0,1,0], 'MarkerSize', 4*opt.size_factor);
    end
    plot(p-offset, SSi_min, '.', 'LineWidth', 1, 'MarkerSize', 13*opt.size_factor, 'Color', 0.9*[1 0 0]);
    plot(p-offset, SSi_max, '.', 'LineWidth', 1, 'MarkerSize', 13*opt.size_factor, 'Color', 0.8*[0 1 0]);
    plot(p-offset, sum(eSi, 3), 'o', 'LineWidth', 0.5, 'MarkerSize', 5*opt.size_factor);
end
if opt.show_expected_final
    plot(p+offset, sum(eSf, 3), 'x', 'LineWidth', 0.5, 'MarkerSize', 5*opt.size_factor);
end

hold off;

if nidx == 1
    title(sprintf('Stored Energy for Storage Unit %d (Gen %d) @ Bus %d', idx, g, b), 'FontSize', 18*opt.size_factor);
else
    if nidx == ns && all(idx == (1:ns)')
        txt = 'All';
    else
        txt = 'Selected';
    end
    title(sprintf('Stored Energy for %s Storage Units', txt), 'FontSize', 18*opt.size_factor);
end
ylabel('Energy, MWh', 'FontSize', 16*opt.size_factor);
xlabel(my_xlabel, 'FontSize', 16*opt.size_factor);
%legend('contract', wlabels{1}, wlabels{2}, wlabels{3}, wlabels{4}, 'Up Res', 'Dn Res', 'Lim', 'Bind', 'Location', 'EastOutside');
set(gca, 'FontSize', 12*opt.size_factor);
h = gcf;
set(h, 'PaperOrientation', 'landscape');
%set(h, 'PaperPosition', [0.25 0.25 10.5 8]);
set(h, 'PaperPosition', [0 0 11 8.5]);
if opt.saveit || opt.saveall
    if nidx == 1
        txt = sprintf('%d', idx);
    elseif nidx == ns && all(idx == (1:ns)')
        txt = 'all';
    else
        txt = 'selected';
    end
    if isempty(strfind(opt.savename, '%s'))
        pdf_name = opt.savename;
    else
        pdf_name = sprintf(opt.savename, txt);
    end
    print('-dpdf', fullfile(opt.savepath, pdf_name));
end
if nargout
    hh = h;
end
