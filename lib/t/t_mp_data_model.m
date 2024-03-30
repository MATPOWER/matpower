function obj = t_mp_data_model(quiet)
% t_mp_data_model - Tests for mp.data_model.

%   MATPOWER
%   Copyright (c) 2020-2024, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

casefile = 't_case_ext';
mpc = loadcase(casefile);
dmc = mp.dm_converter_mpc2().build();
dm0 = mp.data_model().build(mpc, dmc);
nb = size(mpc.bus, 1);
ID2i = zeros(2, 1);     %% initialize as col vector
ID2i(mpc.bus(:, BUS_I)) = (1:nb);

tests = {
    {'filename', casefile},
    {'mpc', mpc},
    {'dm', dm0}
};
nt = length(tests);

t_begin(65*nt + 2, quiet);

for k = 1:nt
    if isa(tests{k}{2}, 'mp.data_model')
        t = sprintf('%s.copy() : ', tests{k}{1});
        dm = tests{k}{2}.copy();
        tests{k}{2}.source = [];
    else
        t = sprintf('mp.data_model().build(%s, dmc) : ', tests{k}{1});
        dm = mp.data_model().build(tests{k}{2}, dmc);
    end
    t_ok(iscell(dm.element_classes), [t 'iscell(dm.element_classes)']);
    t_is(length(dm.element_classes), 5, 12, [t 'length(dm.element_classes)']);
    t_ok(isa(dm.elements, 'mp.mapped_array'), [t 'isa(dm.elements, ''mp.mapped_array'')']);
    t_is(length(dm.elements), 5, 12, [t 'length(dm.elements)']);
    t_ok(isa(dm.elements, 'mp.mapped_array'), [t 'isa(dm.elements, ''mp.mapped_array'')']);
    t_ok(dm.elements.has_name('bus'), [t 'dm.elements.has_name(''bus'')']);
    t_ok(dm.elements.has_name('gen'), [t 'dm.elements.has_name(''gen'')']);
    t_ok(dm.elements.has_name('load'), [t 'dm.elements.has_name(''load'')']);
    t_ok(dm.elements.has_name('branch'), [t 'dm.elements.has_name(''branch'')']);
    t_ok(dm.elements.has_name('shunt'), [t 'dm.elements.has_name(''shunt'')']);
    t_is(dm.elements.name2idx('bus'), 1, 12, [t 'dm.elements.name2idx(''bus'')']);
    t_is(dm.elements.name2idx('gen'), 2, 12, [t 'dm.elements.name2idx(''gen'')']);
    t_is(dm.elements.name2idx('load'), 3, 12, [t 'dm.elements.name2idx(''load'')']);
    t_is(dm.elements.name2idx('branch'), 4, 12, [t 'dm.elements.name2idx(''branch'')']);
    t_is(dm.elements.name2idx('shunt'), 5, 12, [t 'dm.elements.name2idx(''shunt'')']);

    bus = dm.elements.bus;
    t_ok(isa(bus, 'mp.dme_bus'), [t 'bus class']);
    t_ok(isa(bus, 'mp.dm_element'), [t 'bus isa mp.dm_element']);
    t_str_match(bus.name, 'bus', [t 'bus.name']);
    t_is(bus.nr, 10, 12, [t 'bus.nr']);
    t_is(bus.n, 9, 12, [t 'bus.n']);
    t_is(bus.ID, mpc.bus(:, BUS_I), 12, [t 'bus.ID']);
    t_is(bus.ID2i, ID2i, 12, [t 'bus.ID2i']);
    t_is(bus.tab.status, [1;1;1;1;1;0;1;1;1;1], 12, [t 'bus.tab.status']);
    t_is(bus.on, [1;2;3;4;5;7;8;9;10], 12, [t 'bus.on']);
    t_is(bus.off, 6, 12, [t 'bus.off']);

    gen = dm.elements.gen;
    t_ok(isa(gen, 'mp.dme_gen'), [t 'gen class']);
    t_ok(isa(gen, 'mp.dm_element'), [t 'gen isa mp.dm_element']);
    t_str_match(gen.name, 'gen', [t 'gen.name']);
    t_is(gen.nr, 4, 12, [t 'gen.nr']);
    t_is(gen.n, 3, 12, [t 'gen.n']);
    t_is(gen.ID, [1:gen.nr]', 12, [t 'gen.ID']);
    t_is(gen.ID2i, [1:gen.nr]', 12, [t 'gen.ID2i']);
    t_is(gen.tab.status, [1;1;0;1], 12, [t 'gen.tab.status']);
    t_is(gen.on, [1;2;4], 12, [t 'gen.on']);
    t_is(gen.off, 3, 12, [t 'gen.off']);

    ld = dm.elements.load;
    t_ok(isa(ld, 'mp.dme_load'), [t 'load class']);
    t_ok(isa(ld, 'mp.dm_element'), [t 'ld isa mp.dm_element']);
    t_str_match(ld.name, 'load', [t 'ld.name']);
    t_is(ld.nr, 3, 12, [t 'ld.nr']);
    t_is(ld.n, 3, 12, [t 'ld.n']);
    t_is(ld.ID, [1:ld.nr]', 12, [t 'ld.ID']);
    t_is(ld.ID2i, [1:ld.nr]', 12, [t 'ld.ID2i']);
    t_is(ld.tab.status, [1;1;1], 12, [t 'ld.tab.status']);
    t_is(ld.on, [1;2;3], 12, [t 'ld.on']);
    t_ok(isempty(ld.off), [t 'ld.off']);

    branch = dm.elements.branch;
    t_ok(isa(branch, 'mp.dme_branch'), [t 'branch class']);
    t_ok(isa(branch, 'mp.dm_element'), [t 'branch isa mp.dm_element']);
    t_str_match(branch.name, 'branch', [t 'branch.name']);
    t_is(branch.nr, 10, 12, [t 'branch.nr']);
    t_is(branch.n, 9, 12, [t 'branch.n']);
    t_is(branch.ID, [1:branch.nr]', 12, [t 'branch.ID']);
    t_is(branch.ID2i, [1:branch.nr]', 12, [t 'branch.ID2i']);
    t_is(branch.tab.status, [1;1;1;1;1;1;0;1;1;1], 12, [t 'branch.tab.status']);
    t_is(branch.on, [1;2;3;4;5;6;8;9;10], 12, [t 'branch.on']);
    t_is(branch.off, 7, 12, [t 'branch.off']);

    shunt = dm.elements.shunt;
    t_ok(isa(shunt, 'mp.dme_shunt'), [t 'shunt class']);
    t_ok(isa(shunt, 'mp.dm_element'), [t 'shunt isa mp.dm_element']);
    t_str_match(shunt.name, 'shunt', [t 'shunt.name']);
    t_is(shunt.nr, 2, 12, [t 'shunt.nr']);
    t_is(shunt.n, 2, 12, [t 'shunt.n']);
    t_is(shunt.ID, [1:shunt.nr]', 12, [t 'shunt.ID']);
    t_is(shunt.ID2i, [1:shunt.nr]', 12, [t 'shunt.ID2i']);
    t_is(shunt.tab.status, [1;1], 12, [t 'shunt.tab.status']);
    t_is(shunt.on, [1;2], 12, [t 'shunt.on']);
    t_ok(isempty(shunt.off), [t 'shunt.off']);
end

t = 'modify_element_classes : ';
dm = mp.data_model;
e = {@mp.dme_bus, @mp.dme_gen, @mp.dme_load, @mp.dme_branch, @mp.dme_shunt};
t_ok(isequal(dm.element_classes, e), [t 'before']);
dm.modify_element_classes({'mp.dme_shunt', @mp.dme_gizmo, {@mp.dme_gen, 'mp.dme_gen'}});
e = {@mp.dme_bus, @mp.dme_gen, @mp.dme_load, @mp.dme_branch, @mp.dme_gizmo};
t_ok(isequal(dm.element_classes, e), [t 'after']);

t_end;
