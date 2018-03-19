function hh = plot_uc_data(uc1, uc2, optin)
%PLOT_UC_DATA   Plot generator commitment summary
%
%   PLOT_UC_DATA(UC1, UC2, OPT)
%   H = PLOT_UC_DATA(UC1, ...)
%
%   Inputs:
%       UC1     (optional) matrix of 1's and 0's indicating commitment status for
%               first commitment schedule (red)
%               each row corresponds to a generator, each column to a period
%       UC2     (optional) commitment statuses for second commitment schedule
%               (gray)
%       OPT     options struct with the following (all optional) fields, where
%               default values are shown in parenthesis:
%           'title'         ('Unit Commitment - %s') title for the plot, where
%                           %s is an optional placeholder for the subtitle
%           'subtitle'      ({'First', 'Second', 'Both'}) cell array of labels
%                           to use for legend (if two schedules are provided)
%                           and to replace a placeholder in the title (based
%                           on which commitment schedule(s) is(are) provided in
%                           UC1 and/or UC2); can also be a simple string, in
%                           which case no legend will be displayed even if both
%                           UC1 and UC2 are supplied
%           'xlabel'        ('Period') label for horizontal axis
%           'ylabel'        (<empty>) label for vertical axis
%           'rowlabels'     ({'1', '2', '3', ...) labels for rows
%           'saveit'        (false) flag to indicate whether to create PDF file
%           'savepath'      ('') path to directory to save files in
%           'savename'      ('uc-%s.pdf') name of PDF file
%                           %s is optional placeholder for commitment schedule
%                           index ('1', '2', or 'both')
%           'size_factor'   (1) to scale font/marker sizes in case you want to
%                           do sub-plots or something
%
%   Returns handle to current figure window.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% default args
if nargin < 3
    optin = [];
    if nargin < 2
        uc2 = [];
    end
end

%% data and dimensions
if ~isempty(uc1)
    [m, n] = size(uc1);
    if isempty(uc2);
        uc2 = zeros(m, n);
        t = 1;              %% index of subtitle
        txt = '1';          %% for filename
    else
        if size(uc2) ~= [m n]
            error('plot_uc: size mismatch');
        end
        t = 3;              %% index of subtitle
        txt = 'both';       %% for filename
    end
elseif ~isempty(uc2)
    [m, n] = size(uc2);
    uc1 = zeros(m, n);
    t = 2;                  %% index of subtitle
    txt = '2';              %% for filename
end

%% default options
opt = struct( ...
    'title', 'Unit Commitment - %s', ...
    'subtitle', {{'First', 'Second', 'Both'}}, ...
    'xlabel', 'Period', ...
    'ylabel', '', ...
    'rowlabels', [], ...
    'saveit', false, ...
    'savepath', '', ...
    'savename', 'uc-%s.pdf', ...
    'size_factor', 1 );

%% override default options
fields = fieldnames(opt);
for k = 1:length(fields)
    if isfield(optin, fields{k})
        opt.(fields{k}) = optin.(fields{k});
    end
end
if isempty(opt.rowlabels)
    opt.rowlabels = {m:-1:1};
end

%% set subtitle
if iscell(opt.subtitle) && length(opt.subtitle) == 3
    subt = opt.subtitle{t};
elseif ischar(opt.subtitle)
    subt = opt.subtitle;
else
    warning('plot_uc_data: OPT.subtitle must be a string or 3 element cell array of strings');
    subt = '';
end

c1 = [1 0.2 0.2] * 0.8;     %% 1 committed
c2 = [1 1 1] * 0.6;         %% 2 committed
cb = [0.5 0.28 0.28];       %% both committed
cn = [1 1 1];               %% neither committed

h1 = []; h2 = []; hb = [];
for i = 1:m
    for j = 1:n
        y = m-[i-1; i; i; i-1];
        x = [j-1; j-1; j; j];
        if uc1(i, j)                %% 1 committed
            h1 = patch(x, y, c1);
        end
        if uc2(i, j)                %% 2 committed
            h2 = patch(x, y, c2);
        end
        if uc1(i, j) && uc2(i, j)   %% both committed
            hb = patch(x, y, cb);
        end
        if ~uc1(i, j) && ~uc2(i, j) %% neither committed
            hn = patch(x, y, cn);
        end
    end
end

%% legend
if strcmp(txt, 'both') && iscell(opt.subtitle) && length(opt.subtitle) == 3
    legend([h1 h2 hb], opt.subtitle, 'Location', [0.86 0.03 0.1 0.05], 'FontSize', 12*opt.size_factor);
else
    legend('off')
end

%% pretty it up
title(sprintf(opt.title, subt), 'FontSize', 18*opt.size_factor);
axis([0 n 0 m]);
h = gca;
h.XTick = [0.5:1:n-0.5];
h.YTick = [0.5:1:m-0.5];
h.XTickLabel = {1:n};
h.YTickLabel = {m:-1:1};
h.TickLength = [0 0];
h.YAxisLocation = 'right';
h.YTickLabel = opt.rowlabels(m:-1:1);
if ~isempty(opt.xlabel)
    xlabel(opt.xlabel, 'FontSize', 16*opt.size_factor);
end
if ~isempty(opt.ylabel)
    ylabel(opt.ylabel, 'FontSize', 16*opt.size_factor);
end
set(h, 'FontSize', 12*opt.size_factor);

h = gcf;
set(h, 'PaperOrientation', 'landscape');
set(h, 'PaperPosition', [0 0 11 8.5]);
if opt.saveit
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
