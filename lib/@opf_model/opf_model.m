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
%     om = struct( ...
%         'var', struct( ...
%                 'idx', struct( ...
%                         'i1', es, ...
%                         'iN', es, ...
%                         'N', es ), ...
%                 'N', 0, ...
%                 'NS', 0, ...
%                 'order', {}, ...
%                 'data', struct( ...
%                         'v0', es, ...
%                         'vl', es, ...
%                         'vu', es )  ), ...
%         'nln', struct( ...
%                 'idx', struct( ...
%                         'i1', es, ...
%                         'iN', es, ...
%                         'N', es ), ...
%                 'N', 0, ...
%                 'NS', 0, ...
%                 'order', {} ), ...
%         'lin', struct( ...
%                 'idx', struct( ...
%                         'i1', es, ...
%                         'iN', es, ...
%                         'N', es ), ...
%                 'N', 0, ...
%                 'NS', 0, ...
%                 'order', {}, ...
%                 'data', struct( ...
%                         'A', es, ...
%                         'l', es, ...
%                         'u', es )   ), ...
%         'mpc', es   );
    
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
    
    om.mpc = es;

    om = class(om, 'opf_model');
elseif isa(mpc,'opf_model') 
    om = mpc;
else 
    om = opf_model;
    om.mpc = mpc;
end

return;


% om
%     .nvars          total number of variables
%     .nvarsets       number of var sets
%     .nnlcons        total number of non-linear constraints
%     .nnlconsets     number of non-linear constraint sets
%     .nlncons        total number of linear constraints
%     .nlnconsets     number of linear constraint sets
%     .varsets        struct array of var sets
%         .i1         starting index of var set
%         .iN         ending index of var set
%         .N          number of variables in var set
%         .name       name of var set
%         .v0         initial value vector
%         .vl         lower bound vector
%         .vu         upper bound vector
%     .nlconsets      struct array of non-linear constraint sets
%         .i1
%         .iN
%         .N
%         .name
%     .lnconsets      struct array of linear constraint sets
%         .i1
%         .iN
%         .N
%         .name
%         .A
%         .l
%         .u
%     .mpc
%         .baseMVA
%         .bus
%         .branch
%         .gen
%         .gencost
%         .areas
%         .A	(if present, must have l, u)
%         .l
%         .u
%         .N	(if present, must have fparm, H, Cw)
%         .fparm
%         .H
%         .Cw
%     .names
%         .vars
%             .(name)     
%         .nlcons
%         .lncons
