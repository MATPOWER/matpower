function [baseMVA, bus, gen, branch, area, gencost, info] = loadcase(casefile)
%LOADCASE   Load .m or .mat case files or data struct in MATPOWER format
%
%   [baseMVA, bus, gen, branch, area, gencost ] = loadcase(casefile)
%
%   Here casefile is either a struct containing the fields baseMVA, bus,
%   gen, branch, area, gencost, or a string containing the name of the file.
%   If casefile contains the extension '.mat' or '.m', then the explicit file
%   is searched. If casefile containts no extension, then CASELOAD looks for
%   a '.mat' file first, then for a '.m' file.  If the file does not exist
%   or doesn't define all matrices, the routine aborts with an appropriate
%   error message.  Alternatively, it can be called with the syntax:
%
%   [baseMVA, bus, gen, branch, area, gencost, info] = loadcase(casefile)
%
%   In this case, the function will not abort, but info will contain an exit
%   code as follows:
%
%       0:  all variables successfully defined
%       1:  input argument is not a string or struct
%       2:  specified extension-less file name does not exist in search path
%       3:  specified .MAT file does not exist in search path
%       4:  specified .M file does not exist in search path
%       5:  specified file fails to define all matrices

%   MATPOWER
%   $Id$
%   by Carlos Murillo & Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% initialize as empty matrices in case of error
baseMVA = []; bus = []; gen = []; branch = []; area = []; gencost = [];

if isstruct(casefile)				%% first param is a struct
	if nargout < 7
		if ~isfield(casefile,'baseMVA') | ~isfield(casefile,'bus') | ...
				~isfield(casefile,'gen') | ~isfield(casefile,'branch') | ...
				~isfield(casefile,'area') | ~isfield(casefile,'gencost') 
			info = 5;
			return;
		end
	end
	baseMVA	= casefile.baseMVA;
	bus		= casefile.bus;
	gen		= casefile.gen;
	branch	= casefile.branch;
	area	= casefile.area;
	gencost	= casefile.gencost;
	info	= 0;
	return;
else								%% OK, better be a string, else give up
	if ~isstr(casefile)
		if nargout < 7
			error('loadcase: input arg should be a string containing a filename');
	  	else
			info = 1;
			return;
	  	end
	end
end

%% check for explicit extension
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

%% get data from file
if exist('rootname') ~= 1			%% no explicit extension
  if exist([casefile '.mat']) == 2
     load(casefile);
  elseif exist([casefile '.m']) == 2
    [baseMVA, bus, gen, branch, area, gencost] = feval(casefile);
  else
    if nargout < 7
      error('loadcase: specified case not in MATLAB''s search path');
    else
      info = 2;
      return;
    end
  end
else								%% explicit extension given
  if strcmp(extension,'.mat')
    if exist([rootname '.mat']) == 2
      load(rootname) ;
    else
      if nargout < 7
        error('loadcase: specified MAT file does not exist');
      else
        info = 3;
        return;
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
        return;
      end
    end
  end
end

if ~( exist('baseMVA') == 1 & exist('bus') == 1 & exist('branch') == 1 ...
		& exist('gen') == 1 & exist('area') == 1 & exist('gencost') == 1)
  	if nargout < 7
    	error('loadcase: one or more of the data matrices is undefined');
  	else
    	info = 5;
  	end
else
  	info = 0;
end

return;
