function om = display_soln(om, set_type, name, idx)
% display_soln - Display solution values.
% ::
%
%   OM.DISPLAY_SOLN()
%   OM.DISPLAY_SOLN(SET_TYPE)
%   OM.DISPLAY_SOLN(SET_TYPE, NAME)
%   OM.DISPLAY_SOLN(SET_TYPE, NAME, IDX)
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
if nargin < 4
    idx = [];
    if nargin < 3
        name = [];
        if nargin < 2
            set_type = 'all';
        end
    end
end

fid = 1;

mu_thresh = 1e-7;
num_inf = 1e10;

%% print header
if om.is_solved()
    s = om.soln;
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

        if om.(st).N
            if isempty(name)            %% all indices for set type
                idxs = (1:om_st.N);
            elseif isempty(idx)         %% all indices for set type & name
                idxs = [];
                for k = 1:length(om_st.order)
                    if strcmp(om_st.order(k).name, name)
                        i1 = om_st.idx.i1.(name)(om_st.order(k).idx{:});
                        iN = om_st.idx.iN.(name)(om_st.order(k).idx{:});
                        idxs = [idxs (i1:iN)];
                    end
                end
            else                        %% indices for set type, name, idx
                idxs = (om_st.idx.i1.(name)(idx{:}):om_st.idx.iN.(name)(idx{:}));
            end

            %% set type headers
            fprintf(fid, '=====  %s  =====\n', om.set_types.(st));
            %% common columns
            hdr1 = {'  idx    description                ', ...
                    '------- ----------------------------' };
            %% type specific columns
            switch st
            case {'var', 'lin'}
                hdr2 = {'   mu_lb     lb       val      ub      mu_ub', ...
                        ' -------- -------- -------- -------- --------' };
                switch st
                case 'var'
                    [v0, vl, vu] = om.params_var();
                    v = s.x;
                    mu_l = s.lambda.lower;
                    mu_u = s.lambda.upper;
                case 'lin'
                    [A, vl, vu] = om.params_lin_constraint();
                    v = A * s.x;
                    mu_l = s.lambda.mu_l;
                    mu_u = s.lambda.mu_u;
                end
            case 'nle'
                hdr2 = {'    val    lambda', ...
                        ' -------- --------' };
                lam = s.lambda.eqnonlin;
                v = om.eval_nln_constraint(s.x, 1);
            case 'nli'
                hdr2 = {'    val      ub      mu_ub', ...
                        ' -------- -------- --------' };
                mu_u = s.lambda.ineqnonlin;
                v = om.eval_nln_constraint(s.x, 0);
            case 'qdc'
                hdr2 = {'   cost  =  quad    linear  constant  average', ...
                        ' -------- -------- -------- -------- --------' };
                c = [];
                c_k = [];
                c_lin = [];
                c_avg = [];
                vv = om.get_idx();
                for k = 1:length(om_st.order)
                    n = om_st.order(k).name;
                    i = om_st.order(k).idx;
                    [QQ, cc, kk, vs] = om.params_quad_cost(n, i);
                    xx = om.varsets_x(s.x, vs, 'vector');
                    c_total = om.eval_quad_cost(s.x, n, i);
                    len = length(c_total);
                    if len == 1
                        c_constant = kk;
                        c_linear = cc' * xx;
                        if abs(sum(xx)) > 1e-9
                            c_average = c_total / sum(xx);
                        else
                            c_average = NaN;
                        end
                    else
                        c_linear = cc .* xx;
                        c_average = c_total ./ xx;
                        c_average(isinf(c_average)) = NaN;
                        if isscalar(kk)
                            c_constant = ones(len, 1)*kk/len;
                        else
                            c_constant = kk;
                        end
                    end
                    if sum(sum(QQ)) == 0 && sum(kk) == 0 && len == length(cc)
                        c_average = cc;
                    end
                    c = [c; c_total];
                    c_lin = [c_lin; c_linear];
                    c_k = [c_k; c_constant];
                    c_avg = [c_avg; c_average];
                end
                c_quad = c - c_lin - c_k;
            case 'nlc'
                c = [];
                for k = 1:length(om_st.order)
                    n = om_st.order(k).name;
                    i = om_st.order(k).idx;
                    c_total = om.eval_nln_cost(s.x, n, i);
                    c = [c; c_total];
                end
                hdr2 = {'   cost', ...
                        ' --------' };
            otherwise
                hdr2 = {'', ...
                        '' };
            end

            %% print headers
            for h = 1:length(hdr1)
                fprintf(fid, '%s\n', [hdr1{h} hdr2{h}]);
            end

            %% print data
            none = '- ';
            % none = '';
            for k = 1:length(idxs)
                ii = sprintf('%d', idxs(k));
                fmt = sprintf('%%-%ds', length(ii)+ceil((7-length(ii))/2));
                fprintf(fid, '%7s %-28s', sprintf(fmt, ii), om.describe_idx(st, idxs(k)));

                switch st
                case {'var', 'lin'}
                    if abs(mu_l(idxs(k))) < mu_thresh
                        mu_lb = sprintf(none);
                    else
                        mu_lb = sprintf_num(8, mu_l(idxs(k)));
                    end
                    if abs(mu_u(idxs(k))) < mu_thresh
                        mu_ub = sprintf(none);
                    else
                        mu_ub = sprintf_num(8, mu_u(idxs(k)));
                    end
                    if vl(idxs(k)) < -num_inf
                        lb = sprintf(none);
                    else
                        lb = sprintf_num(8, vl(idxs(k)));
                    end
                    if vu(idxs(k)) > num_inf
                        ub = sprintf(none);
                    else
                        ub = sprintf_num(8, vu(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s%9s%9s\n', ...
                        mu_lb, lb, sprintf_num(8, v(idxs(k))), ub, mu_ub);
                case 'nle'
                    if abs(lam(idxs(k))) < mu_thresh
                        mu_ub = sprintf(none);
                    else
                        mu_ub = sprintf_num(8, lam(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s\n', sprintf_num(8, v(idxs(k))), ...
                        mu_ub);
                case 'nli'
                    if abs(mu_u(idxs(k))) < mu_thresh
                        mu_ub = sprintf(none);
                    else
                        mu_ub = sprintf_num(8, mu_u(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s\n', sprintf_num(8, v(idxs(k))), ...
                        '0', mu_ub);
                case 'qdc'
                    if isnan(c_avg(idxs(k)))
                        cc_avg = '- ';
                    else
                        cc_avg = sprintf_num(8, c_avg(idxs(k)));
                    end
                    fprintf(fid, '%9s%9s%9s%9s%9s\n', ...
                        sprintf_num(8, c(idxs(k))), ...
                        sprintf_num(8, c_quad(idxs(k))), ...
                        sprintf_num(8, c_lin(idxs(k))), ...
                        sprintf_num(8, c_k(idxs(k))), ...
                        cc_avg);
                case 'nlc'
                    fprintf(fid, '%9s\n', sprintf_num(8, c(idxs(k))));
                otherwise
                    fprintf(fid, '\n');
                end
            end     %% loop over idxs

            %% section footer
            switch st
            case {'var', 'lin'}
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Min', ...
                    sprintf_num(8, min(mu_l)), sprintf_num(8, min(vl)), ...
                    sprintf_num(8, min(v)), ...
                    sprintf_num(8, min(vu)), sprintf_num(8, min(mu_u)));
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', 'Max', ...
                    sprintf_num(8, max(mu_l)), sprintf_num(8, max(vl)), ...
                    sprintf_num(8, max(v)), ...
                    sprintf_num(8, max(vu)), sprintf_num(8, max(mu_u)));
                fprintf(fid, '\n');
            case 'nle'
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s\n', '', 'Min', ...
                    sprintf_num(8, min(v)), ...
                    sprintf_num(8, min(lam)));
                fprintf(fid, '%7s %-28s%9s%9s\n', '', 'Max', ...
                    sprintf_num(8, max(v)), ...
                    sprintf_num(8, max(lam)));
                fprintf(fid, '\n');
            case 'nli'
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s\n', '', 'Max', ...
                    sprintf_num(8, max(v)), '0', ...
                    sprintf_num(8, max(mu_u)));
                fprintf(fid, '\n');
            case 'qdc'
                fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', ...
                    'Sum of Displayed Costs', ...
                    sprintf_num(8, sum(c(idxs))), ...
                    sprintf_num(8, sum(c_quad(idxs))), ...
                    sprintf_num(8, sum(c_lin(idxs))), ...
                    sprintf_num(8, sum(c_k(idxs))), '');
                fprintf(fid, '\n');
            case 'nlc'
                if length(c) > 1
                    fprintf(fid, '%s\n', [hdr1{2} hdr2{2}]);
                    fprintf(fid, '%7s %-28s%9s%9s%9s%9s%9s\n', '', ...
                        'Sum of Displayed Costs', ...
                        sprintf_num(8, sum(c(idxs))));
                end
                fprintf(fid, '\n');
            otherwise
                fprintf(fid, '\n');
            end
        end
    end             %% loop over set types
else
    fprintf(fid, 'Not a solved model.\n');
end

function str = sprintf_num(width, val)
val = full(val);
fmt = sprintf('%%%dg', width);
str = sprintf(fmt, val);
if length(str) > width
    fmt = sprintf('%%%d.*g', width);
    for p = width-2:-1:0
        str = sprintf(fmt, p, val);
        if length(str) <= width
            break;
        end
    end
end
assert(length(str) <= width)
