function hh = plot_uc(md, varargin)
%PLOT_UC   Plot generator commitment summary
%
%   PLOT_UC(MD)
%   PLOT_UC(MD, IDX)
%   PLOT_UC(MD, IDX, '<option1_name>', '<option1_value', ...)
%   PLOT_UC(MD, IDX, OPT)
%   PLOT_UC(MD, IDX, OPT, '<option1_name>', '<option1_value', ...)
%   PLOT_UC(MD, MD2)
%   PLOT_UC(MD, MD2, IDX, )
%   PLOT_UC(MD, MD2, IDX, '<option1_name>', '<option1_value', ...)
%   PLOT_UC(MD, MD2, IDX, OPT)
%   PLOT_UC(MD, MD2, IDX, OPT, '<option1_name>', '<option1_value', ...)
%   H = PLOT_UC(MD, ...)
%
%   IDX is a vector of gen indices, if empty, it includes all generators in
%   MD. Options can be specified as an OPT struct or as individual pairs of
%   'name' and 'value' arguments. The possible options include the following,
%   where the default is shown in parenthesis:
%       'title'         ('Unit Commitment - %s') title for the plot, where
%                       %s is an optional placeholder for the subtitle
%       'subtitle'      ({'First', 'Second', 'Both'}) cell array of labels
%                       to use for legend (if two schedules are provided)
%                       and to replace a placeholder in the title (based
%                       on which commitment schedule(s) is(are) provided in
%                       UC1 and/or UC2); can also be a simple string, in
%                       which case no legend will be displayed even if both
%                       UC1 and UC2 are supplied
%       'xlabel'        ('Period') label for horizontal axis
%       'ylabel'        (<empty>) label for vertical axis
%       'rowlabels'     ({'1', '2', '3', ...) labels for rows (top to bottom)
%       'saveit'        (false) flag to indicate whether to create PDF file
%       'saveall'       (false) flag to indicate whether to create a single
%                       PDF file or, if both MD and MD2 are supplied, three
%                       PFF files, one for MD, one for MD2 and one for both.
%       'savepath'      ('') path to directory to save files in
%       'savename'      ('uc-%s.pdf') name of PDF file
%                       %s is optional placeholder for storage unit index
%       'size_factor'   (1) to scale font/marker sizes in case you want to
%                       do sub-plots or something
%
%   Returns handle to current figure window.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
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
nt = md.idx.nt;
ng = size(md.mpc.gen, 1);
nj_max = max(md.idx.nj);
nc_max = max(max(md.idx.nc));

%% default optional args
idx = [];
md2 = [];
k = 1;
if nargin > 1                   %% 2nd arg is ...
    if isstruct(varargin{k})    %% ... md2
        md2 = varargin{k};
    else                        %% ... idx
        idx = varargin{k};
    end
    k = k + 1;
end
if isempty(idx) && nargin > 2   %% 3rd arg is ...
    if isnumeric(varargin{k})
        idx = varargin{k};      %% ... idx
        k = k + 1;
    end
end

%% default idx
if isempty(idx)
    idx = find(any(md.UC.CommitKey == 0, 2) | any(md.UC.CommitKey == 1, 2));
end
nidx = length(idx);
if nidx > 1 && size(idx, 1) == 1
    idx = idx';     %% convert row vector to column vector
end
b = md.mpc.gen(idx, GEN_BUS);

%% default options
opt = struct( ...
    'title', 'Unit Commitment - %s', ...
    'subtitle', {{'First', 'Second', 'Both'}}, ...
    'xlabel', 'Period', ...
    'ylabel', '', ...
    'rowlabels', [], ...
    'saveit', false, ...
    'saveall', false, ...
    'savepath', '', ...
    'savename', 'uc-%s.pdf', ...
    'size_factor', 1 );

%% process options
vargs = varargin(k:end); %% remaining arguments
if mod(length(vargs), 2)    %% odd number of options, first must be OPT struct
    if ~isstruct(vargs{1})
        error('plot_uc: Single OPT argument must be a struct');
    end
    myopt = vargs{1};
    k = 2;
else                        %% even number of options
    myopt = struct;
    k = 1;
end
while k < length(vargs)
    opt_name = vargs{k};
    opt_val  = vargs{k+1};
    if ~isfield(opt, opt_name)
        error('plot_uc: ''%s'' is not a valid option name', opt_name);
    end
    myopt.(opt_name) = opt_val;
    k = k + 2;
end
fields = fieldnames(myopt);
for f = 1:length(fields)
    opt.(fields{f}) = myopt.(fields{f});
end

%% check CommitKey match
uc1 = md.UC.CommitSched(idx, :);
if ~isempty(md2)
    uck2 = md2.UC.CommitKey(idx, :);
    if any(any( uck2 ~= md.UC.CommitKey(idx, :) & (uck2 == 2 | uck2 == 0) ))
        error('plot_uc: CommitKey fields in MD and MD2 do not match');
    end
    uc2 = md2.UC.CommitSched(idx, :);
else
    uc2 = [];
end

%% generator labels
m = size(uc1, 1);
opt.rowlabels = cell(m, 1);
if isfield(md.mpc, 'genfuel')
    for i = 1:m
        opt.rowlabels{i} = sprintf('Gen %d @ Bus %d, %s', idx(i), b(i), md.mpc.genfuel{idx(i)});
    end
else
    for i = 1:m
        opt.rowlabels{i} = sprintf('Gen %d @ Bus %d', idx(i), b(i));
    end
end

%% do the plots
if opt.saveall
    opt.saveit = 1;
end
if opt.saveall || isempty(uc2)
    clf;
    h = plot_uc_data(uc1, [], opt);
end
if ~isempty(uc2)
    if opt.saveall || isempty(uc1)
        clf;
        h = plot_uc_data([], uc2, opt);
    end
    if opt.saveall || ~isempty(uc1)
        clf;
        h = plot_uc_data(uc1, uc2, opt);
    end
end

if nargout
    hh = h;
end
