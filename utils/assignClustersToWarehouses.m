function [assignIdx, locker2WarehouseAssignment] = assignClustersToWarehouses(clusterCenters, ...
    warehouseCoords, idx, varargin)
%ASSIGNCLUSTERSTOWAREHOUSES Greedy one-to-one mapping from clusters to warehouses.
%C
%   [ASSIGNIDX, TOTALCOST, D] = ASSIGNCLUSTERSTOWAREHOUSES(C, W)
%   assigns each cluster (rows of C) to a distinct warehouse (rows of W)
%   using the greedy rule:
%     - take clusters in order (1..K by default),
%     - assign each cluster to the nearest warehouse that is not yet assigned.
%
%   Inputs:
%     clusterCenters  - Kc-by-d matrix (cluster centroids)
%     warehouseCoords - Kw-by-d matrix (warehouse coordinates)
%     idx             - nLockers-by-1 vector, idx(i)=j if locker i is
%                           assingned to cluster j
%
%   Optional name/value pairs:
%     'ClusterOrder' - permutation of 1:Kc specifying the order to process
%                      clusters (default 1:Kc)
%     'Metric'       - distance metric for pdist2 (default 'euclidean')
%
%   Outputs:
%     assignIdx                     - Kc-by-1 vector where assignIdx(i) = j means cluster i -> warehouse j
%     locker2WarehouseAssignment    - nLockers-by-1 vector, v(i)=j if locker i is
%                                       assigned to warehouse j
%
%   Notes:
%     - The function requires Kw >= Kc (enough warehouses for a unique assignment).
%     - If you want a different processing order, pass 'ClusterOrder'.
%
%   Example:
%     C = [0 0; 1 1; 5 5];        % 3 cluster centroids
%     W = [0 1; 2 2; 6 6; 10 10]; % 4 warehouses (>= clusters)
%     [idx,cost] = assignClustersToWarehouses_greedy(C, W);
%     disp(idx); disp(cost);
%

% parse inputs
p = inputParser;
addRequired(p, 'clusterCenters', @(x) isnumeric(x) && ndims(x)==2);
addRequired(p, 'warehouseCoords', @(x) isnumeric(x) && ndims(x)==2);
addParameter(p, 'ClusterOrder', [], @(x) isempty(x) || (isvector(x) && all(x==floor(x))));
addParameter(p, 'Metric', 'euclidean', @(s) ischar(s) || isstring(s));
parse(p, clusterCenters, warehouseCoords, varargin{:});

C = p.Results.clusterCenters;
W = p.Results.warehouseCoords;
order = p.Results.ClusterOrder;
metric = char(p.Results.Metric);

Kc = size(C,1);
Kw = size(W,1);

if Kw < Kc
    error('Number of warehouses (%d) must be >= number of clusters (%d) for unique greedy assignment.', Kw, Kc);
end

% default order
if isempty(order)
    order = 1:Kc;
else
    if numel(order) ~= Kc || ~all(sort(order(:)') == 1:Kc)
        error('ClusterOrder must be a permutation of 1:Kc.');
    end
end

% distance matrix (clusters x warehouses)
distMatrix = pdist2(C, W, metric); % Kc x Kw

assignIdx = zeros(Kc,1);        % cluster -> warehouse
assignedWarehouses = false(Kw,1);

% Greedy assignment in specified order
for ii = 1:Kc
    ci = order(ii);
    drow = distMatrix(ci, :);
    drow(assignedWarehouses) = inf;   % ignore already assigned warehouses
    [~, wh] = min(drow);
    assignIdx(ci) = wh;
    assignedWarehouses(wh) = true;
end

nLockers = length(idx);
locker2WarehouseAssignment = zeros(nLockers, 1);
for i=1:nLockers
    clusterIdx = idx(i);
    warehouseIdx = assignIdx(clusterIdx);
    locker2WarehouseAssignment(i) = warehouseIdx;
end

end
