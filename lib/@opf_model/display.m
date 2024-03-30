function display(om)
% display - Displays the object.
% ::
%
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
% See also opt_model.

%   MATPOWER
%   Copyright (c) 2008-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% call parent with added set type
display@opt_model(om, {'cost'});

fprintf('  mpc = ');
if ~isempty(fieldnames(om.mpc))
    fprintf('\n');
end
if have_feature('octave')
    fprintf('    <scalar struct>\n');
else
    display(om.mpc);
end
