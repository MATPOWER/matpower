function [groups, unvisited] = connected_components(C, groups, unvisited)
%CONNECTED_COMPONENTS Returns the connected components of a graph
%   [GROUPS, ISOLATED] = CONNECTED_COMPONENTS(C)
%
%   Returns the connected components of a directed graph, specified by
%   a node-branch incidence matrix C, where C(I, J) = -1 if node J is
%   connected to the beginning of branch I, 1 if it is connected to
%   the end of branch I, and zero otherwise. The return value GROUPS
%   is a cell array of vectors of the node indices for each component.
%   The second return value ISOLATED is a vector of indices of isolated
%   nodes that have no connecting branches.

%   Can also be called with a current GROUP and list of as-yet
%   UNVISITED nodes:
%   [GROUPS, UNVISITED] = CONNECTED_COMPONENTS(C, GROUPS, UNVISITED)
%   Internally, this is used recursively.

%   MATPOWER
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% initialize groups and unvisited list
nn = size(C, 2);        %% number of nodes
Ct = C';
visited = zeros(1, nn); %% start with nothing visited
if nargin < 2
    groups = {};
    unvisited = 1:nn;
    isolated = find(sum(abs(C), 1) == 0);   %% find isolated nodes
    unvisited(isolated) = [];   %% remove isolated nodes from unvisited
end

%% traverse graph, starting with first unvisited node
cn = unvisited(1);      %% set current node to first unvisited node
visited(cn) = 1;        %% mark current node as visited
%% use FIFO queue for breadth-first, LIFO stack for depth-first
qs = zeros(1, nn);      %% initialize queue/stack
f = 1; N = 0;           %% queue starts at position f, has N elements
% t = 0;                  %% top of stack at position t

%% push current node to queue/stack
qs(mod(f+N-1, nn)+1) = cn; N = N + 1;       %% FIFO / BFS
% t = t + 1; qs(t) = cn;                      %% LIFO / DFS

while N                 %% FIFO / BFS
% while t                 %% LIFO / DFS
    %% pop next node
    cn = qs(f); N = N - 1; f = mod(f, nn) + 1;  %% FIFO / BFS
%     cn = qs(t); t = t - 1;                      %% LIFO / DFS

    %% find all nodes connected to current node
    %% (finding rows of column-indexed Ct, rather than cols of row-indexed C,
    %%  because row-indexing a sparse matrix is sloooowww, results in ~30x
    %%  speed up on ~60k bus network)
    [jj, junk] = find(Ct(:, C(:, cn) ~= 0));    %% non-zeros in rows connected to cn
%    [jj, ~] = find(Ct(:, C(:, cn) ~= 0));   %% non-zeros in rows connected to cn
%    [~, jj] = find(C(C(:, cn) ~= 0, :));    %% non-zeros in rows connected to cn
    cnn = jj(visited(jj) == 0); %% indices of non-visited cols (may contain dups)

    %% mark them as visited and queue them
    for k = 1:length(cnn)
        if visited(cnn(k)) == 0         %% if not yet visited
            visited(cnn(k)) = 1;        %% mark as visited
            %% push it to queue/stack
            N = N + 1; qs(mod(f+N-2, nn)+1) = cnn(k);   %% FIFO / BFS
%             t = t + 1; qs(t) = cnn(k);                  %% LIFO / DFS
        end
    end
end

%% add visited nodes to group
group = find(visited);
groups{end+1} = group;

%% update unvisited
v = ones(nn, 1);
v(unvisited) = 0;
v(group)     = 1;
unvisited    = find(v == 0);

%% recurse to traverse remaining components
if isempty(unvisited)       %% sort groups by cardinality
    l = cellfun('length', groups);
    [junk, i] = sort(l, 2, 'descend');
%     [~, i] = sort(l, 2, 'descend');
    groups = groups(i);
else                        %% recurse
    [groups, unvisited] = connected_components(C, groups, unvisited);
end

%% prepare to return isolated nodes
if nargin < 2
    unvisited = isolated;
end
