function mpc = order_radial(mpc)
%ORDER_RADIAL  Performs oriented ordering to buses and branches.
%
%   mpc = order_radial(mpc)
%
%   orders the branches by the the principle of oriented ordering:
%   indicies of sending nodes are smaller than the indicies of the
%   receiving nodes. The branch index is equal to the index of their
%   receiving node. The ordering is taken from:
%   D. Rajicic, R. Ackovski and R. Taleski, "Voltage correction power flow,"
%   IEEE Transactions on Power Delivery, vol. 9, no. 2, pp. 1056-1062, Apr 1994.
%
%   See also RADIAL_PF.

%% define named indices into bus, gen, branch matrices
define_constants;
%% Initialize
slack = mpc.bus(mpc.bus(:,BUS_TYPE) == REF, 1);
[f, t] = deal(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS));
nl = size(mpc.branch,1);
branch_number = (1:nl)';
mpc.branch_order = [];
mpc.loop = [];
mpc.bus_order = slack;
%% Bus and branch ordering
iter = 1;
while ~isempty(f) && iter <= nl
    % For each of the "from" buses that are in bus_order,
    % add the corresponding "to" buses to it. If both "from" and "to" are
    % in bus_order, the branch is forming loop. Add corresponding branch
    % numbers to branch_order or to loop.
    mf = ismember(f,mpc.bus_order);
    mt = ismember(t,mpc.bus_order);
    is_loop = mf & mt;
    % Add branch to loop and delete it from the list
    if any(is_loop)
        mpc.loop = [mpc.loop; branch_number(is_loop)];
        mf(is_loop) = [];
%         mt(is_loop) = [];
        f(is_loop) = [];
        t(is_loop) = [];
        branch_number(is_loop) = [];
    end
    % Add unique buses to bus_order
    if any(mf)
        u = unique(t(mf)); % unique buses
        [junk,i] = intersect(t.*mf,u); % first occurence of unique buses in t
        mpc.bus_order = [mpc.bus_order; t(i)];
        % Add branch to branch_order and delete it from the list
        mpc.branch_order = [mpc.branch_order; branch_number(i)];
        mf(i) = [];
%         mt(i) = [];
        f(i) = [];
        t(i) = [];
        branch_number(i) = [];
    end
    % Add any remaining branch to loop and delete it from the list
    if any(mf)
        mpc.loop = [mpc.loop; branch_number(mf)];
        f(mf) = [];
        t(mf) = [];
        branch_number(mf) = [];
    end
    % For each of the "to" buses that are in bus_order,
    % add the corresponding "from" buses to it. If both "from" and "to" are
    % in bus_order, the branch is forming loop. Add corresponding branch
    % numbers to branch_order or to loop.
    mf = ismember(f,mpc.bus_order);
    mt = ismember(t,mpc.bus_order);
    is_loop = mf & mt;
    % Add branch to loop and delete it from the list
    if any(is_loop)
        mpc.loop = [mpc.loop; branch_number(is_loop)];
        mt(is_loop) = [];
%         mf(is_loop) = [];
        f(is_loop) = [];
        t(is_loop) = [];
        branch_number(is_loop) = [];
    end
    % Add unique buses to bus_order
    if any(mt)
        u = unique(f(mt)); % unique buses
        [junk,i] = intersect(f.*mt,u); % first occurence of unique buses in f
        mpc.bus_order = [mpc.bus_order; f(i)];
        % Add branch to branch_order and delete it from the list
        mpc.branch_order = [mpc.branch_order; branch_number(i)];
%         mf(i) = [];
        mt(i) = [];
        f(i) = [];
        t(i) = [];
        branch_number(i) = [];
    end
    % Add any remaining branch to loop and delete it from the list
    if any(mt)
        mpc.loop = [mpc.loop; branch_number(mt)];
        f(mt) = [];
        t(mt) = [];
        branch_number(mt) = [];
    end
    iter = iter + 1;
end
if ~isempty(f)
    mpc.not_connected = branch_number;
else
    mpc.not_connected = [];
end
%% Reorder bus, branch and gen.
if isempty(mpc.loop)
    % Permute rows in branch
    mpc.branch = mpc.branch(mpc.branch_order,:);
    % Make an "inverse" vector out of bus_order
    mpc.bus_order_inv = sparse(mpc.bus_order,ones(nl+1,1),1:nl+1);
    % Swap indicies in "from" and "to" buses using bus_order_inv
    [f, t] = deal(mpc.branch(:,F_BUS),mpc.branch(:,T_BUS));
    f = mpc.bus_order_inv(f);
    t = mpc.bus_order_inv(t);
    % Reverse branch orientation of "from" is biger than "to"
    mpc.br_reverse = f > t;
    [f(mpc.br_reverse), t(mpc.br_reverse)] = deal(t(mpc.br_reverse), f(mpc.br_reverse));
    % Put new "from" and "to" indicies in branch
    mpc.branch(:,[F_BUS T_BUS]) = [f t];
    % Make an "inverse" vector out of branch_order
    mpc.branch_order_inv = sparse(mpc.branch_order,ones(nl,1),1:nl);
    % Permute rows in bus and replace bus indicies
    mpc.bus = mpc.bus(mpc.bus_order,:);
    mpc.bus(:,BUS_I) = mpc.bus_order_inv(mpc.bus(:,BUS_I));
    % Replace bus indicies in gen
    mpc.gen(:,GEN_BUS) = mpc.bus_order_inv(mpc.gen(:,GEN_BUS));
end