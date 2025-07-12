function m = milp_ex1
% milp_ex1 - Example of mixed-integer linear program (MILP) optimization.
%
% Example of using MP-Opt-Model to build and solve an optimization model,
% in this case, a mixed-integer, linear programming (MILP) problem.
%
% **Multi-Plant Production and Distribution Example Problem**
% 
% A company operates two production plants and needs to satisfy customer
% demand for two products, Y and Z.
%
%   - Plants can be opened or closed.
%   - Each plant incurs a fixed cost if it is used.
%   - Each plant has a capacity limit (maximum total units it can ship).
%   - Products must be shipped from plants to customers.
%   - Products are produced as needed for shipping, with no extra inventory.
%   - Delivery costs depend on the plant, customer, and product and include
%     production and shipping.
%   - All customer demand must be fully satisfied (no shortages allowed).
%   - Shipping from a plant is only allowed if that plant is open.
%
% Goals
%
% The goal is to decide:
%
%   - Which plants to open.
%   - How much of each product to ship from each open plant to each customer.
%
% To minimize total cost:
%
%   - Sum of plant fixed costs (if used)
%   - Plus delivery costs for all products shipped

%%-----  Problem Data  -----
PlantCapacity = [ 60; 60 ];     % Plants 1, 2

% Customer Demand for each Product: total = 56 units
CustomerDemandY = [ 10; 15; 5 ];    % Customers 1, 2, 3
CustomerDemandZ = [  8; 12; 6 ];    % Customers 1, 2, 3

% Delivery costs for Product Y and Z
DeliveryCostY = [
%   Plant 1     Plant 2
    4           7;          % Customer 1
    5           6;          % Customer 2
    3           8           % Customer 3
];
DeliveryCostZ = [
%   Plant 1     Plant 2
    10          4;          % Customer 1
     9          5;          % Customer 2
    12          3           % Customer 3
];

% Plant fixed cost for 3 scenarios:
%   1 : plant 2 relatively more expensive, shut down
%   2 : plant 1 relatively more expensive, shut down
%   3 : both plants operating
PlantFixedCost = [
%   Scenario 1  Scenario 2  Scenario 3
    100         100          50;    % Plant 1
    200         100         120     % Plant 2
];

%%-----  Build Model  -----
Scenario = 1;   % build scenario 1 initially

% create the model object
mm = mp.opt_model();

%% Decision variables
% plant status
% add set of binary variable of dimension 2, named 'u'
mm.var.add('u', 2, [], 0, 1, 'B');

% product delivery
% prepare to add 2 sets of 'y' and 'z' variables, indexed by plant
mm.var.init_indexed_name('y', {2});
mm.var.init_indexed_name('z', {2});

% add them, each with dimension 3, one for each customer
% bounded below by 0, no upper bound
for p = 1:2
    mm.var.add('y', {p}, 3, [], 0);     % Product Y, Plant p to Customers 1-3
    mm.var.add('z', {p}, 3, [], 0);     % Product Z, Plant p to Customers 1-3
end

% display model variables
% show starting/ending indices, dimensions for each set
% mm.var

%% Costs
% plant fixed costs
% add fixed costs, providing coefficients for 'u' variables only
mm.qdc.add(mm.var, 'fixed', [], PlantFixedCost(:, Scenario), [], {'u'});

% delivery costs
% prepare to add 2 sets of delivery costs per product, indexed by plant
mm.qdc.init_indexed_name('delivery_y', {2});
mm.qdc.init_indexed_name('delivery_z', {2});

% add costs for plant p, for each product, corresponding to variable sets
% 'y{p}' and 'z{p}'; each is 3 x 1 corresponding to the 3 customers
for p = 1:2
    vs_y = struct('name', 'y', 'idx', {{p}});    % variable set for 'y{p}'
    vs_z = struct('name', 'z', 'idx', {{p}});    % variable set for 'z{p}'
    mm.qdc.add(mm.var, 'delivery_y', {p}, [], DeliveryCostY(:, p), [], vs_y);
    mm.qdc.add(mm.var, 'delivery_z', {p}, [], DeliveryCostZ(:, p), [], vs_z);
end

% display model costs
% show starting/ending indices, dimensions for each set
% mm.qdc

%% Constraints
% plant capacity constraints
% define based on full set of variables: u y{1} y{2} z{1} z{2}
Au = -spdiags(PlantCapacity, 0, 2, 2);
Ayz = sparse([1 1 1 1 1 1 0 0 0 0 0 0;   % sum for plant 1
              0 0 0 0 0 0 1 1 1 1 1 1]); % sum for plant 2
ub = 0;    % constraint upper bound (no lower bound)
mm.lin.add(mm.var, 'capacity', [Au Ayz], [], ub);

% demand satisfaction constraints
% demand for product Y: based on variable sets: y{1} y{2} only
vs_y = struct('name', 'y', 'idx', {{1}, {2}});
Ad = [speye(3) speye(3)];
mm.lin.add(mm.var, 'demand_y', Ad, CustomerDemandY, CustomerDemandY, vs_y);

% demand for product Z: based on variable sets: z{1} z{2} only
vs_z = struct('name', 'z', 'idx', {{1}, {2}});
mm.lin.add(mm.var, 'demand_z', Ad, CustomerDemandZ, CustomerDemandZ, vs_z);

% display model constraints
% show starting/ending indices, dimensions for each set
% mm.lin

%%-----  Solve Model  -----
opt = struct('verbose', 2, 'skip_prices', true);
mm.solve(opt);

% display solution
mm.display_soln();

%%-----  Three Scenarios  -----
for Scenario = 1:3
    mm.qdc.set_params(mm.var, 'fixed', 'c', PlantFixedCost(:, Scenario));
    mm.solve();

    y1 = sum(mm.var.get_soln(mm.soln, 'y', {1}));   % product Y from plant 1
    y2 = sum(mm.var.get_soln(mm.soln, 'y', {2}));   % product Y from plant 2
    z1 = sum(mm.var.get_soln(mm.soln, 'z', {1}));   % product Z from plant 1
    z2 = sum(mm.var.get_soln(mm.soln, 'z', {2}));   % product Z from plant 2

    subplot(1, 3, Scenario);
    bar([y1 z1; y2 z2], 'stacked');
    ax = axis();
    axis([ax(1:3) 60]);
    title(sprintf('Scenario %d, Total Cost = %g', Scenario, mm.soln.f));
    legend({'Product Y', 'Product Z'});
    xlabel('Plants');
    ylabel('Production (units)');

    % mm.var.display_soln(mm.soln);
end

if nargout
    m = mm;
end
