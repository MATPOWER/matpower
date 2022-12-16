function tab_class = mp_table_class()
%MP_TABLE_CLASS  Returns handle to constructor for TABLE or MP_TABLE
%
%   Returns a handle to TABLE constructor, if it is available, otherwise to
%   MP_TABLE constructor.
%
%   See also TABLE, MP_TABLE.

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
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
