function om = display_soln(om, varargin)
% display_soln - Display solution values.
% ::
%
%   OM.DISPLAY_SOLN()
%   OM.DISPLAY_SOLN(SET_TYPE)
%   OM.DISPLAY_SOLN(SET_TYPE, NAME)
%   OM.DISPLAY_SOLN(SET_TYPE, NAME, IDX)
%   OM.DISPLAY_SOLN(FID)
%   OM.DISPLAY_SOLN(FID, SET_TYPE)
%   OM.DISPLAY_SOLN(FID, SET_TYPE, NAME)
%   OM.DISPLAY_SOLN(FID, SET_TYPE, NAME, IDX)
%
%   Displays the model solution, including values, bounds and shadow
%   prices for variables and linear constraints, values and shadow
%   prices for nonlinear constraints, and individual cost components.
%
%   Results are displayed for each SET_TYPE or specified SET_TYPE and
%   for each named/indexed set or a specified NAME/IDX.
%
%   Inputs:
%       SET_TYPE - one of the following, specifying the type of set:
%           'var' - variables
%           'lin' - linear constraints
%           'nle' - nonlinear equality constraints
%           'nli' - nonlinear inequality constraints
%           'nlc' - nonlinear costs
%           'qdc' - quadratic costs
%         or
%           a cell array of one or more of the above
%         or
%           '' or 'all' - indicating to display all
%       NAME - (optional) char array specifying the name of the set
%       IDX  - (optional) cell array specifying the indices of the set
%
%   Examples:
%       om.display_soln('var');
%       om.display_soln({'nle', 'nli'});
%       om.display_soln('var', 'P');
%       om.display_soln('lin', 'lin_con_1');
%       om.display_soln('nle', 'nle_con_b', {2,3});
%
% See also get_soln, parse_soln.

%   MP-Opt-Model
%   Copyright (c) 2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MP-Opt-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/mp-opt-model for more info.

%% input arg handling
if nargin < 2 || ischar(varargin{1})
    fid = 1;
    args = varargin;
else
    fid = varargin{1};
    args = varargin(2:end);
end
nargs = length(args);

if nargs < 3
    idx = [];
    if nargs < 2
        name = [];
        if nargs < 1
            set_type = 'all';
        end
    end
end

%% print header
if om.is_solved()
    if strcmp(set_type, 'all')
        set_types = fieldnames(om.set_types);   %% all set types
    elseif ~iscell(set_type)
        set_types = {set_type}; %% make set_type a cell array of char arrays
    else
        set_types = set_type;
    end

    for ss = 1:length(set_types)
        st = set_types{ss};
        om_st = om.(st);

        switch st
        case 'var'
            om_st.display_soln(om.soln, fid, args{2:end});
        case 'nle'
            om_st.display_soln(om.var, om.soln, 1, fid, args{2:end});
        case 'nli'
            om_st.display_soln(om.var, om.soln, 0, fid, args{2:end});
        otherwise
            om_st.display_soln(om.var, om.soln, fid, args{2:end});
        end
    end             %% loop over set types
else
    fprintf(fid, 'Not a solved model.\n');
end
