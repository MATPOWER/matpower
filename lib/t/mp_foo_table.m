classdef mp_foo_table < mp_table_subclass

%   MATPOWER
%   Copyright (c) 2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function obj = mp_foo_table(varargin)
            obj@mp_table_subclass(varargin{:});
        end
    end     %% methods
end         %% classdef
