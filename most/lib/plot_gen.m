function hh = plot_gen(md, idx, varargin)
%PLOT_GEN   Plot generator results
%
%   PLOT_GEN(MD)
%   PLOT_GEN(MD, IDX)
%   PLOT_GEN(MD, IDX, '<option1_name>', '<option1_value', ...)
%   PLOT_GEN(MD, IDX, OPT)
%   PLOT_GEN(MD, IDX, OPT, '<option1_name>', '<option1_value', ...)
%   H = PLOT_GEN(MD, ...)
%
%   IDX is the gen index. If IDX is a vector, it sums them first,
%   if empty, it includes all generators in MD. Options can be
%   specified as an OPT struct or as individual pairs of 'name' and 'value'
%   arguments. The possible options include the following, where the default
%   is shown in parenthesis:
%       'saveit'        (false) flag to indicate whether to create PDF file
%       'saveall'       (false) flag to indicate whether to create individual
%                       PDF files for each element of IDX, as well as
%                       aggregate, when IDX is a vector (or empty)
%       'savepath'      ('') path to directory to save files in
%       'savename'      ('gen-%s.pdf') name of PDF file
%                       %s is optional placeholder for storage unit index
%       'size_factor'   (1) to scale font/marker sizes in case you want to
%                       do sub-plots or something
%       'show_Pc'                   (true) Pc, the energy contract
%       'show_variable_Pmax'        (true) non-constant Pmax
%       'show_limits'               (true) max Pmax and min Pmin
%       'show_contingencies'        (true) contingency dispatches
%       'show_reserves'             (false) contingency & ramping reserves
%       'show_grid'                 (true) vertical lines to divide periods
%
%   Returns handle to current figure window.

%   TO DO: Include separate plot/subplot of reserves, ramping reserves.
%          Include options for enabling/disabling specific portions.
%          Do initial invisible plot to get v axis parameters.
%          Make sure it makes sense for dispatchable loads (should
%          indicate dispatchable load in title).

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

%% gather data
my_xlabel = 'Period';
nt = md.idx.nt;
ng = size(md.mpc.gen, 1);
nj_max = max(md.idx.nj);
nc_max = max(max(md.idx.nc));

%% input args
if nargin < 2
    idx = [];
end
if isempty(idx)
    idx = (1:ng)';
end
nidx = length(idx);
if nidx > 1 && size(idx, 1) == 1
    idx = idx';     %% convert row vector to column vector
end
b = md.mpc.gen(idx, GEN_BUS);

%% default options
opt = struct( ...
    'saveit', false, ...
    'saveall', false, ...
    'savepath', '', ...
    'savename', 'gen-%s.pdf', ...
    'tol', 1e-3, ...
    'size_factor', 1, ...
    'show_grid', true, ...
    'show_Pc', true, ...
    'show_variable_Pmax', 1, ...
    'show_limits', 1, ...
    'show_contingencies', 1, ...
    'show_reserves', 0 );

%% process options
if mod(length(varargin), 2) %% odd number of options, first must be OPT struct
    if ~isstruct(varargin{1})
        error('plot_gen: Single OPT argument must be a struct');
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
        error('plot_gen: ''%s'' is not a valid option name', opt_name);
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
        error('plot_gen: ''savename'' must include a ''%s'' placeholder when ''saveall'' option is true.');
    end
    for i = 1:nidx
        plot_gen(md, idx(i), opt, 'saveit', true);
    end
end

%% initialize data structures to be plotted
if isfield(md.mpc, 'genfuel')
    gf = md.mpc.genfuel{idx};
else
    gf = '';
end
p = (1:nt)';
Pg = zeros(nt, nj_max);
Pgk = NaN(nt, nj_max, nc_max);
maxPmax = -Inf(nt, 1);
minPmin =  Inf(nt, 1);
Pmax = NaN(nt, nj_max);
ncPmax = NaN(nt, nj_max);
ncPg = NaN(nt, nj_max);

%% Pmax constant over j?
constantPmax = 1;
constantPmin = 1;
for t = 1:nt
    maxPmax(t) = sum(md.flow(t,1,1).mpc.gen(idx, PMAX), 1);
    minPmin(t) = sum(md.flow(t,1,1).mpc.gen(idx, PMIN), 1);
    for j = 1:md.idx.nj(t)
        if sum(md.flow(t,j,1).mpc.gen(idx, PMAX), 1) ~= maxPmax(t)
            constantPmax = 0;
            if maxPmax(t) < sum(md.flow(t,j,1).mpc.gen(idx, PMAX), 1);
                maxPmax(t) = sum(md.flow(t,j,1).mpc.gen(idx, PMAX), 1);
            end
        end
        if sum(md.flow(t,j,1).mpc.gen(idx, PMIN), 1) ~= minPmin(t)
            constantPmin = 0;
            if minPmin(t) > sum(md.flow(t,j,1).mpc.gen(idx, PMIN), 1);
                minPmin(t) = sum(md.flow(t,j,1).mpc.gen(idx, PMIN), 1);
            end
        end
        Pg(t,j) = sum(md.flow(t,j,1).mpc.gen(idx, PG), 1);
        for k = 1:md.idx.nc(t,j)
            Pgk(t,j,k) = sum(md.flow(t,j,k+1).mpc.gen(idx, PG), 1);
        end
        Pmax(t,j) = sum(md.flow(t,j,1).mpc.gen(idx, PMAX), 1);
        if abs(Pmax(t,j) - Pg(t,j)) > opt.tol
            ncPmax(t,j) = Pmax(t,j);
            ncPg(t,j) = Pg(t,j);
        end
    end
end
ePg = sum(md.results.ExpectedDispatch(idx, :), 1);
Pc = sum(md.results.Pc(idx, :), 1)';
Rpp = sum(md.results.Rpp(idx, :), 1)';
Rpm = sum(md.results.Rpm(idx, :), 1)';
Rrp = sum(md.results.Rrp(idx, :), 1)';
Rrm = sum(md.results.Rrm(idx, :), 1)';
Gmax = Pc + Rpp;
Gmin = Pc - Rpm;

%% do plots
%% figure out axis limits
plot(p, maxPmax+0.001, 'LineStyle', 'none')
hold on
plot(p, minPmin-0.001, 'LineStyle', 'none')

%% draw central path (patch)
maxPg = max(Pg, [], 2)+0.001;
minPg = min(Pg, [], 2)-0.001;
patch([p;p(end:-1:1); p(1)], [minPg; maxPg(end:-1:1); minPg(1)], 0.95*[1 1 1], 'LineStyle', 'none')

%% adjust vertical axis if necessary
v = axis;
if v(3) > 0             %% make sure zero is included on vertical axis
    v(3) = 0;
elseif v(4) < 0
    v(4) = 0;
end
if opt.show_reserves    %% expand for reserves if necessary
    if v(4) < max([max(Rpp), max(Rpm), max(Rrp), max(Rrm)])
        v(4) = max([max(Rpp), max(Rpm), max(Rrp), max(Rrm)]) + 0.001;
    end
end

%% draw grid
if opt.show_grid
    m = v(3) - 100*(v(4) - v(3));
    M = v(4) + 100*(v(4) - v(3));
    for k = 0:nt
        line([k; k], [m; M], 'LineWidth', 0.25, 'Color', 0.9*[1 1 1]);
    end
end
axis(v);

%% dispatches
if opt.show_reserves
    plot(p, Rpp, '-.', 'Color', 0.8*[0 1 0], 'LineWidth', 1);
    plot(p, Rpm, '-.', 'Color', 0.9*[1 0 0], 'LineWidth', 1);
    plot(p, [NaN; Rrp], ':', 'Color', 0.8*[0 1 0], 'LineWidth', 1);
    plot(p, [NaN; Rrm], ':', 'Color', 0.9*[1 0 0], 'LineWidth', 1);
end
plot(p, Pg, '-', 'Color', 0.8*[1 1 1], 'LineWidth', 1);
if opt.show_Pc
    plot(p, Pc, '--', 'Color', 0.5*[1 1 1], 'LineWidth', 2);
end
plot(p, Gmax, '--', 'Color', 0.8*[0 1 0], 'LineWidth', 1);
plot(p, Gmin, '--', 'Color', 0.9*[1 0 0], 'LineWidth', 1);
plot(p, ePg', 'Color', [0 0 1], 'LineWidth', 2);
if nc_max && opt.show_contingencies
    plot(p, reshape(Pgk, nt, nj_max*nc_max), 'x');
end
plot(p, Pg, 'o', 'MarkerSize', 6*opt.size_factor);

%% limits
if opt.show_limits
    plot(p, maxPmax, ':k', 'LineWidth', 1);
    plot(p, minPmin, ':k', 'LineWidth', 1);
end
if ~constantPmax && opt.show_variable_Pmax
%     plot(p, ncPmax, 'v', 'Color', 0.8*[0 1 0], 'LineWidth', 1, 'MarkerSize', 6*opt.size_factor);
    plot(p, Pmax, 'v', 'MarkerSize', 6*opt.size_factor);
%     plot(p, ncPg, '^', 'Color', 0.9*[1 0 0], 'LineWidth', 1, 'MarkerSize', 6*opt.size_factor);
%     plot(p, ncPg, '^', 'LineWidth', 1, 'MarkerSize', 6*opt.size_factor);
end

hold off;

if nidx == 1
    if ~isempty(md.Storage.UnitIdx)
        s = find(idx == md.Storage.UnitIdx);
    else
        s = 0;
    end
    if s
        title(sprintf('Real Power Output for Storage Unit %d (Gen %d) @ Bus %d', s, idx, b), 'FontSize', 18*opt.size_factor);
    else
        tt = sprintf('Real Power Output for Gen %d @ Bus %d', idx, b);
        if ~isempty(gf)     %% add fuel type if available
            tt = sprintf('%s (%s)', tt, gf);
        end
        title(tt, 'FontSize', 18*opt.size_factor);
    end
else
    if nidx == ng && all(idx == (1:ng)')
        txt = 'All';
    else
        txt = 'Selected';
    end
    title(sprintf('Real Power Output for %s Generators', txt), 'FontSize', 18*opt.size_factor);
end
ylabel('Real Power, MW', 'FontSize', 16*opt.size_factor);
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
    elseif nidx == ng && all(idx == (1:ng)')
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
