function om = opf_model(mpc)
%OPF_MODEL  Constructor for OPF model class
%
% om
%   .var
%       .idx
%           .i1
%           .iN
%           .N
%       .N
%       .NS
%       .data
%           .v0
%           .vl
%           .vu
%       .order
%   .nln
%       .idx
%           .i1
%           .iN
%           .N
%       .N
%       .NS
%       .order
%   .lin
%       .idx
%           .i1
%           .iN
%           .N
%       .N
%       .NS
%       .data
%           .A
%           .l
%           .u
%           .vs
%       .order
%   .cost
%       .idx
%           .i1
%           .iN
%           .N
%       .N
%       .NS
%       .data
%           .N
%           .H
%           .Cw
%           .dd
%           .rr
%           .kk
%           .mm
%           .vs
%       .order
%   .mpc
%       .baseMVA
%       .bus
%       .branch
%       .gen
%       .gencost
%       .areas
%       .A	(if present, must have l, u)
%       .l
%       .u
%       .N	(if present, must have fparm, H, Cw)
%       .fparm
%       .H
%       .Cw

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% es = struct();    %% doesn't work in Matlab 6
es = struct('tmp', 0);
es = rmfield(es, 'tmp');
if nargin == 0
    om.var.idx.i1 = es;
    om.var.idx.iN = es;
    om.var.idx.N = es;
    om.var.N = 0;
    om.var.NS = 0;
    om.var.order = {};
    om.var.data.v0 = es;
    om.var.data.vl = es;
    om.var.data.vu = es;

    om.nln.idx.i1 = es;
    om.nln.idx.iN = es;
    om.nln.idx.N = es;
    om.nln.N = 0;
    om.nln.NS = 0;
    om.nln.order = {};

    om.lin.idx.i1 = es;
    om.lin.idx.iN = es;
    om.lin.idx.N = es;
    om.lin.N = 0;
    om.lin.NS = 0;
    om.lin.order = {};
    om.lin.data.A = es;
    om.lin.data.l = es;
    om.lin.data.u = es;
    om.lin.data.vs = es;
    
    om.cost.idx.i1 = es;
    om.cost.idx.iN = es;
    om.cost.idx.N = es;
    om.cost.N = 0;
    om.cost.NS = 0;
    om.cost.order = {};
    om.cost.data.N = es;
    om.cost.data.H = es;
    om.cost.data.Cw = es;
    om.cost.data.dd = es;
    om.cost.data.rh = es;
    om.cost.data.kk = es;
    om.cost.data.mm = es;
    om.cost.data.vs = es;
    om.cost.params = es;
    
    om.mpc = es;

    om = class(om, 'opf_model');
elseif isa(mpc,'opf_model') 
    om = mpc;
else 
    om = opf_model;
    om.mpc = mpc;
end

return;
