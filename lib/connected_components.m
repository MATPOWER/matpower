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
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2012 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as MATLAB code and MEX-files) available in a
%   MATLAB(R) or comparable environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%% initialize groups and unvisited list
nn = size(C, 2);        %% number of branchs and nodes
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
queue = cn;             %% enqueue current node
while ~isempty(queue)
    %% dequeue next node
    cn = queue(end);  queue(end) = [];  %% use FIFO queue for breadth-first
%   cn = queue(1);    queue(1)   = [];  %% use LIFO stack for depth-first

    %% find all nodes connected to current node
    [~, jj] = find(C(C(:, cn) ~= 0, :));    %% non-zeros in rows connected to cn
    cnn = jj(visited(jj) == 0); %% indices of non-visited cols (may contain dups)

    %% mark them as visited and queue them
%   visited(cnn) = 1;           %% mark as visited
%   queue = [cnn; queue];       %% enqueue them
    for k = 1:length(cnn)
        if visited(cnn(k)) == 0         %% if not yet visited
            visited(cnn(k)) = 1;        %% mark as visited
            queue = [cnn(k); queue];    %% enqueue it
        end
    end
end

%% add visited nodes to group
group = find(visited);
groups{end+1} = group;
% ng = length(group);
% [m, n] = size(groups);
% if ng > n       %% expand groups before appending new group
%     groups = [ groups zeros(m, ng-n); group ];
% else
%     groups = [ groups; group zeros(1, n-ng) ];
% end

%% update unvisited
v = ones(nn, 1);
v(unvisited) = 0;
v(group)     = 1;
unvisited    = find(v == 0);

%% recurse to traverse remaining components
if ~isempty(unvisited)
    [groups, unvisited] = connected_components(C, groups, unvisited);
end

%% prepare to return isolated nodes
if nargin < 2
    unvisited = isolated;
end
