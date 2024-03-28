function TorF = is_mixed_integer(om)
% is_mixed_integer - Return true if model is mixed integer, false otherwise.
% ::
%
%   TorF = OM.IS_MIXED_INTEGER()
%
%   Outputs:
%       TorF : 1 or 0, indicating whether any of the variables are
%              binary or integer
%
% See also opt_model.

%   MP-Opt-Model
%   Copyright (c) 2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%   To do: Make this a property that gets set true whenever binary or
%          integer vars are added to the problem.

TorF = 0;
if om.getN('var')
    for k = 1:length(om.var.order)
        t = om.var.data.vt.(om.var.order(k).name);
        if iscell(t)
            for j = 1:length(t(:))
                if any(t{j} ~= 'C')
                    TorF = 1;
                    break;
                end
            end
        else
            if any(t ~= 'C')
                TorF = 1;
                break;
            end
        end
    end
end