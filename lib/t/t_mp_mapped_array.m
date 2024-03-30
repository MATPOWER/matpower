function obj = t_mp_mapped_array(quiet)
% t_mp_mapped_array - Tests for mp.mapped_array.

%   MATPOWER
%   Copyright (c) 2021-2024, Power Systems Engineering Research Center (PSERC)
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

t_begin(69, quiet);

%% set up example data
vals = {[1;2;3], {'one','two','three'}, struct('a',10,'b',20), mp.mapped_array()};
names = {'vector', 'cellstr', 'struct', 'object'};
% m0 = cell2struct(num2cell(1:4), names, 2);
N = length(vals);

t = 'A = mp.mapped_array() : ';
A = mp.mapped_array();
t_ok(isa(A, 'mp.mapped_array'), [t 'class']);
t_is(length(A), 0, 12, [t 'length']);

t = 'A = mp.mapped_array(vals) : ';
A = mp.mapped_array(vals);
t_ok(isa(A, 'mp.mapped_array'), [t 'class']);
t_is(length(A), 4, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals), [t 'vals']);
t_ok(isequal(A.p_.names, cell(1, N)), [t 'names']);

t = 'A = mp.mapped_array(vals, names) : ';
A = mp.mapped_array(vals, names);
t_ok(isa(A, 'mp.mapped_array'), [t 'class']);
t_is(length(A), 4, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals), [t 'vals']);
t_ok(isequal(A.p_.names, names), [t 'names']);

t = 'has_name : ';
t_ok(A.has_name('vector'), [t 'vector']);
t_ok(A.has_name('cellstr'), [t 'cellstr']);
t_ok(A.has_name('struct'), [t 'struct']);
t_ok(A.has_name('object'), [t 'object']);
t_ok(~A.has_name('random'), [t 'random']);

t = 'name2idx : ';
t_is(A.name2idx('vector'), 1, 12, [t 'vector']);
t_is(A.name2idx('cellstr'), 2, 12, [t 'cellstr']);
t_is(A.name2idx('struct'), 3, 12, [t 'struct']);
t_is(A.name2idx('object'), 4, 12, [t 'object']);

t = 'subsref : ';
t_ok(isequal(A.cellstr, vals{2}), [t 'b = A.<name>']);
t_ok(isequal(A{2}, vals{2}), [t 'b = A{i}']);
if have_feature('octave')   %% https://savannah.gnu.org/bugs/?60726
                            %% https://savannah.gnu.org/bugs/?48693
    t_skip(1, 'b = {A{ii<vector>}} does not work due to Octave incompatibility');
else
    t_ok(isequal({A{[3,1]}}, {vals{[3,1]}}), [t 'b = {A{idxs}}']);
end
t_ok(isequal(A([4;2]), vals([4;2])), [t 'b = A(idxs)']);

t_ok(isequal(A.struct.a, 10), [t 'b = A.<name>.<name2>']);
t_ok(isequal(A.cellstr{3}, 'three'), [t 'b = A.<name>{i}']);
t_ok(isequal(A.vector(1), 1), [t 'b = A.<name>(i)']);
t_ok(isequal(A{3}.b,  20), [t 'b = A{i}.<name>']);
t_ok(isequal(A{2}{2}, 'two'), [t 'b = A{i}{j}']);
t_ok(isequal(A{1}(3), 3), [t 'b = A{i}(j)']);

A.vector = [3;2;1];
t_ok(isequal(A{1}, [3;2;1]), [t 'A.<name> = b']);
A{3} = struct('x', 1, 'y', 2);
t_ok(isequal(A.struct, struct('x', 1, 'y', 2)), [t 'A{i} = b']);
A([4;2]) = vals([2;4]);
t_ok(isequal(A([2;3;4]), {vals{4},A.struct,vals{2}}), [t 'A(ii<vector>) = b']);

t = 'subsasgn : ';
A.struct.x = 10;
t_ok(isequal(A{3}.x, 10), [t 'A.<name>.<name2> = b']);
A.object{1} = 'uno';
t_ok(isequal(A{4}{1}, 'uno'), [t 'A.<name>{i} = b']);
A.vector([4;1]) = 4;
t_ok(isequal(A{1}(1:4), [4;2;1;4]), [t 'A.<name>(i) = b']);
A{3}.z = 30;
t_ok(isequal(A.struct.z, 30), [t 'A{i}.<name2> = b']);
A{4}{3} = 'tres';
t_ok(isequal(A.object{3}, 'tres'), [t 'A{i}{j} = b']);
A{1}(:) = [2:2:8];
t_ok(isequal(A.vector, [2;4;6;8]), [t 'A{i}(i) = b']);

t = 'delete_elements(i) : ';
A = mp.mapped_array(vals, names);
A.delete_elements(1);
t_is(length(A), 3, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals(2:4)), [t 'vals']);
t_ok(isequal(A.p_.names, names(2:4)), [t 'names']);
m = struct('cellstr', 1, 'struct', 2, 'object', 3);
t_ok(isequal(A.p_.map, m), [t 'map']);

t = 'delete_elements(<name>) : ';
A = mp.mapped_array(vals, names);
A.delete_elements('cellstr');
t_is(length(A), 3, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals([1 3 4])), [t 'vals']);
t_ok(isequal(A.p_.names, names([1 3 4])), [t 'names']);
m = struct('vector', 1, 'struct', 2, 'object', 3);
t_ok(isequal(A.p_.map, m), [t 'map']);

t = 'delete_elements(ii<vector>) : ';
A = mp.mapped_array(vals, names);
A.delete_elements([3 1]);
t_is(length(A), 2, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals([2 4])), [t 'vals']);
t_ok(isequal(A.p_.names, names([2 4])), [t 'names']);
m = struct('cellstr', 1, 'object', 2);
t_ok(isequal(A.p_.map, m), [t 'map']);

t = 'delete_elements({names}) : ';
A = mp.mapped_array(vals, names);
A.delete_elements({'vector', 'object'});
t_is(length(A), 2, 12, [t 'length']);
t_ok(isequal(A.p_.vals, vals([2 3])), [t 'vals']);
t_ok(isequal(A.p_.names, names([2 3])), [t 'names']);
m = struct('cellstr', 1, 'struct', 2);
t_ok(isequal(A.p_.map, m), [t 'map']);

t = 'add_elements({<single>}) : ';
A = mp.mapped_array(vals, names);
A.add_elements({{'hello', 'world'}});
t_is(length(A), 5, 12, [t 'length']);
t_ok(isequal(A{5}, {'hello', 'world'}), [t 'vals']);

t = 'add_elements({<multiple>}) : ';
A.add_elements({'hi', pi});
t_is(length(A), 7, 12, [t 'length']);
t_ok(isequal(A(6:7), {'hi', pi}), [t 'vals']);

t = 'add_elements({<single>}, names) : ';
A.add_elements({struct('x',1/2,'y',1/4)}, {'coord'});
t_is(length(A), 8, 12, [t 'length']);
t_ok(isequal(A{8}, struct('x',1/2,'y',1/4)), [t 'vals']);
t_ok(isequal(A.p_.names{8}, 'coord'), [t 'names']);

t = 'add_elements({<multiple>}, names) : ';
A.add_elements({exp(1), cell(2,2)}, {'e', 'emptycell'});
t_is(length(A), 10, 12, [t 'length']);
t_ok(isequal(A(9:10), {exp(1), cell(2,2)}), [t 'vals']);
t_ok(isequal(A.p_.names(9:10), {'e', 'emptycell'}), [t 'names']);

t = 'add_names';
A.add_names(5, {'hello', 'saludos', 'pi'});
t_ok(isequal(A.p_.names(5:7), {'hello', 'saludos', 'pi'}), [t 'names']);

t = 'B = A.copy() : ';
A = mp.mapped_array(vals, names);
A.object = mp.mapped_array(vals, names);
B = A.copy();
A.delete_elements('vector');
A.object.cellstr{2} = 'changed';
t_is(length(A), 3, 12, [t 'length(A)']);
t_is(length(B), 4, 12, [t 'length(B)']);
t_ok(isequal(A.object.cellstr{2}, 'changed'), [t 'A has modified object field']);
t_ok(isequal(B.object.cellstr{2}, 'two'), [t 'B has original object field']);

if nargout
    obj = A;
end

t_end;
