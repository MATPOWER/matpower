function tab_class = mp_table_class()
% mp_table_class - Returns handle to constructor for :class:`table` or mp_table.
%
% Returns a handle to :class:`table` constructor, if it is available,
% otherwise to mp_table constructor.
%
% ::
%
%   % Works in MATLAB or Octave, which does not support table().
%   table_class = mp_table_class();
%   T = table_class(var1, var2, ...);
%
% See also table, mp_table.

%   MATPOWER
%   Copyright (c) 2022-2023, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if have_feature('table')
    tab_class = @table;
else
    tab_class = @mp_table;
end