function mv = mipsver
%MIPSVER  Prints or returns MIPS version information.
%   V = MIPSVER returns the current MIPS version number. Calling MIPSVER
%   without assigning the return value prints the version and release date of
%   the current installation of MIPS.

%   MIPS
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/minopf/ for more info.

v = struct( 'Name',     'MIPS', ... 
            'Version',  '1.0', ...
            'Release',  '', ...
            'Date',     '08-Mar-2010' );
if nargout > 0
    mv = v.Version;
else
    fprintf('%-22s Version %-9s  %11s\n', v.Name, v.Version, v.Date);
end
