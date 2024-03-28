function display(om, varargin)
% display - Displays the object.
%
% Called when semicolon is omitted at the command-line. Displays the details
% of the variables, constraints, costs included in the model.
%
% See also opt_model.

%   MP-Opt-Model
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% Due to a bug related to inheritance in constructors in
%% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
%% INIT_SET_TYPES() cannot be called directly in the
%% MP_IDX_MANAGER constructor, as desired.
%%
%% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
%%              after object construction, but before object use.
if isempty(om.var)          %% only if not already initialized
    om.init_set_types();
end

if nargin < 2
    more_set_types = {};
else
    more_set_types = varargin{1};
end

%% display details of each set type
set_types = {'var', 'nle', 'nli', 'lin', 'qdc', 'nlc', more_set_types{:}};
fprintf('\n');
for k = 1:length(set_types)
    om.display_set(set_types{k});
end

%% user data
fields = fieldnames(om.userdata);
if ~isempty(fields)
    fprintf('\nUSER DATA\n')
    fprintf('=========\n')
    fprintf('  name                               size       class\n');
    fprintf(' ------------------------------   -----------  --------------------\n');
    for k = 1:length(fields)
        f = om.userdata.(fields{k});
        [m, n] = size(f);
        fprintf('  %-31s %5dx%-5d   %s\n', fields{k}, m, n, class(f));
    end
else
    fprintf('USER DATA                   :  <none>\n');
end
