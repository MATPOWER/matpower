function display(om)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also OPT_MODEL.

%   MP-Opt-Model
%   Copyright (c) 2008-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% display details of each set type
set_types = {'var', 'nle', 'nli', 'lin', 'qdc', 'nlc'};
set_names = struct(...
    'var',  'VARIABLES', ...
    'nle',  'NONLIN EQ CONSTRAINTS', ...
    'nli',  'NONLIN INEQ CONSTRAINTS', ...
    'lin',  'LINEAR CONSTRAINTS', ...
    'qdc',  'QUADRATIC COSTS', ...
    'nlc',  'GEN NONLIN COSTS'  );
fprintf('\n');
for k = 1:length(set_types)
    om.display_set(set_types{k}, set_names.(set_types{k}));
end

%% user data
fprintf('  userdata = ');
if ~isempty(om.userdata)
    fprintf('\n');
end
if have_feature('octave')
    fprintf('    <scalar struct>\n');
else
    display(om.userdata);
end
