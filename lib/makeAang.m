function [Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt)
%

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% options
ignore_ang_lim = mpopt(25);     %% OPF_IGNORE_ANG_LIM

%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% data dimensions
nl = size(branch, 1);   %% number of branches

if ignore_ang_lim
  Aang 	= sparse(0, nb);
  lang 	= [];
  uang 	= [];
  iang 	= [];
else
  iang = find((branch(:, ANGMIN) & branch(:, ANGMIN) > -360) | ...
			  (branch(:, ANGMAX) & branch(:, ANGMAX) < 360));
  iangl = find(branch(iang, ANGMIN));
  iangh = find(branch(iang, ANGMAX));
  nang = length(iang);

  if nang > 0
	ii = [(1:nang)'; (1:nang)'];
	jj = [branch(iang, F_BUS); branch(iang, T_BUS)];
	Aang = sparse(ii, jj, [ones(nang, 1); -ones(nang, 1)], nang, nb);
	uang = Inf * ones(nang,1);
	lang = -uang;
	lang(iangl) = branch(iang(iangl), ANGMIN) * pi/180;
	uang(iangh) = branch(iang(iangh), ANGMAX) * pi/180;
  else
	Aang = sparse(0, nb);
	lang =[];
	uang =[];
  end
end

return;
