classdef (Abstract) form < handle
%MP.FORM  MATPOWER Formulation abstract base class.
%   Each concrete Network Model Element class must inherit, at least
%   indirectly, from both MP.NM_ELEMENT and MP.FORM.
%
%   MP.FORM provides properties and methods related to the specific
%   formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       subclasses provide properties for model parameters
%
%   Methods
%       form_name() - returns string w/name of formulation
%       form_tag() - returns string w/short label for formulation
%       model_params() - cell array of names of model parameters
%       get_params() - 
%       find_form_class() - 
%       superclass_tab() - 

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         mp_form_field = '';
%     end

    methods
        function name = form_name(obj)
            error('mp.form/form_name: must be implemented in subclass');
        end

        function tag = form_tag(obj)
            error('mp.form/form_tag: must be implemented in subclass');
        end

        function params = model_params(obj)
            error('mp.form/model_params: must be implemented in subclass');
        end

        function varargout = get_params(obj, idx, names)
            % [p1, p2, ..., pN] = obj.get_params(idx)
            % pA = obj.get_params(idx, nameA)
            % [pA, pB, ...] = obj.get_params(idx, {nameA, nameB, ...})
            if nargin < 3
                names = obj.model_params();
            end
            if nargin < 2
                idx = [];
            end
            if ~iscell(names)
                names = {names};
            end
            np = obj.nk * obj.np;
            ncols = [1; np; obj.nk * obj.nz];   %% 1 = 1, 2 = np, 3 = nz

            if nargout > length(names)
                namestr = sprintf(' ''%s''', names{:})
                error('mp.form/get_params: return values must correspond to%s', namestr);
            end

            varargout = cell(1, nargout);
            if isempty(idx)     %% all ports
                for k = 1:nargout
                    P = obj.(names{k});
                    if isempty(P)
                        varargout{k} = sparse(np, ncols(obj.param_ncols.(names{k})));
                    else
                        varargout{k} = P;
                    end
                end
            else                %% selected ports
                ni = length(idx);
                for k = 1:nargout
                    P = obj.(names{k});
                    if isempty(P)
                        varargout{k} = sparse(ni, ncols(obj.param_ncols.(names{k})));
                    else
                        varargout{k} = P(idx, :);
                    end
                end
            end
        end

        function form_class = find_form_class(obj)
            if isa(obj, 'mp.form')
                tab = obj.superclass_tab({'mp.form', 'mp.nm_element'});

                %% select among classes that are not mp.nm_element's ...
                j = find(tab.ii(:, 2) == 0);
                %% ... the mp.form class with longest inheritance path
                %% back to mp.form
                [~, k] = max(tab.ii(j, 1));

                assert(~isempty(k));
                form_class = tab.name{j(k)};
            else
                form_class = '<none>';
            end
        end

        function [tab, ii] = superclass_tab(obj, roots, mcls, tab, ii, level)
            %% obj - an mp.form object
            %% roots - cell array of names of root classes to search for
            %% mcls - meta class object (undocumented MATLAB and Octave) for
            %%        the class of interest
            %% tab - struct with fields:
            %%      name - cell array of class names
            %%      ii   - matrix where each row corresponds to entry in name
            %%             field, and column j corresponds to j-th entry in
            %%             roots, value is generations of inheritance required
            %%             to reach that root
            %% ii - row of tab.ii corresponding to mcls on input (before
            %%      traversing parent classes)
            %% level - indent level for this class when printing
            %%         (hard-coded 'verbose' variable not 0)
            verbose = 0;
            if nargin < 6
                level = 0;
                if nargin < 5
                    ii = [];
                    if nargin < 4
                        tab = struct('name', {{}}, 'ii', []);
                        if nargin < 3
                            mcls = meta.class.fromName(class(obj));
                            if nargin < 2
                                roots = {};
                            end
                        end
                    end
                end
            end
            n = length(roots);
            if isempty(ii)
                ii = zeros(1, n);
            end
            done = zeros(1, n);
            i0 = ii;
            for k = 1:n
                if strcmp(mcls.Name, roots{k})
                    done(k) = 1;
                end
            end
            if verbose
                prefix = repmat(' ', 1, 4*level);
                fprintf('%s %s\n', prefix, mcls.Name);
            end
            if have_feature('octave')
                sclist = mcls.SuperClassList;
            else
                sclist = {};
                for k = 1:length(mcls.SuperclassList)
                    sclist{end+1} = mcls.SuperclassList(k);
                end
            end
            for k = 1:length(sclist)
                [tab, iii] = obj.superclass_tab(roots, sclist{k}, tab, i0, level+1);
                for k = 1:n
                    if done(k) == 1
                        ii(k) = 1;
                    elseif iii(k) && i0(k) < iii(k) + 1
                        ii(k) = iii(k) + 1;
                    end
                end
            end
            if verbose
                fmt = repmat(' %d', 1, n);
                fprintf(['%s' fmt '\n'], prefix, ii);
            end
            tab.name{end+1, 1} = mcls.Name;
            tab.ii(end+1, 1:n) = ii;
        end
    end     %% methods
end         %% classdef
