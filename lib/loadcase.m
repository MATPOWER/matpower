function [baseMVA, bus, gen, branch, area, gencost, info] = loadcase(casefile)
%CASELOAD   Load .m or .mat case files in MATPOWER format
%
%   [baseMVA, bus, gen, branch, area, gencost ] = loadcase(casefile)
%
%   Here casefile is a string containing the name of the file.  If casefile
%   contains the extension '.mat' or '.m', then that explicit file is
%   searched. If casefile containts no extension, then CASELOAD looks for
%   a '.mat' file first, then for a '.m' file.  If the file does not exist
%   or doesn't define all matrices, the routine aborts with an appropriate
%   error message.  Alternatively, it can be called with the syntax:
%
%   [baseMVA, bus, gen, branch, area, gencost, info] = loadcase(casefile)
%
%   In this case, the function will not abort, but info will contain an exit
%   code as follows:
%
%       0:  all variables succesfully defined
%       1:  input argument is not a string
%       2:  specified extension-less file name does not exist in search path
%       3:  specified .MAT file does not exist in search path
%       4:  specified .M file does not exist in search path
%       5:  specified file fails to define all matrices

%   MATPOWER
%   $Id$
%   by Carlos Murillo, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

if ~isstr(casefile)
  if nargout < 7
    error('loadcase: input arg should be a string containing a filename');
  else
    info = 1;
    baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];
    return
  end
end

l = length(casefile);
if l > 2
  if strcmp(casefile(l-1:l), '.m')
    rootname = casefile(1:l-2);
    extension = '.m';
  elseif l > 4
    if strcmp(casefile(l-3:l), '.mat')
       rootname = casefile(1:l-4);
       extension = '.mat';
    end
  end
end

if exist('rootname') ~= 1
  if exist([casefile '.mat']) == 2
     load(casefile);
  elseif exist([casefile '.m']) == 2
    [baseMVA, bus, gen, branch, area, gencost] = feval(casefile);
  else
    if nargout < 7
      error('loadcase: specified case not in MATLAB''s search path');
    else
      info = 2;
      baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];
      return
    end
  end
else
  if strcmp(extension,'.mat')
    if exist([rootname '.mat']) == 2
      load(rootname) ;
    else
      if nargout < 7
        error('loadcase: specified MAT file does not exist');
      else
        info = 3;
        baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];
        return
      end
    end
  elseif strcmp(extension,'.m')
    if exist([rootname '.m']) == 2
      [baseMVA, bus, gen, branch, area, gencost] = feval(rootname);
    else
      if nargout < 7
        error('loadcase: specified M file does not exist');
      else
        info = 4;
        baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];
        return
      end
    end
  end
end

if ~( exist('baseMVA') == 1 & exist('bus') == 1 & exist('branch') == 1 ...
     & exist('area') == 1 & exist('gencost') == 1)
  if nargout < 7
    error('loadcase: one or more of the data matrices is undefined');
  else
    info = 1;
    baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];
  end
else
  info = 0;
end

return;
