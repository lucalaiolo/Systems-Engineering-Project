function [idx, centroids] = clusterizeLockers(X, k, varargin)
% Balanced k-means clustering (equal-ish cluster sizes).
%
%   [IDX, CENTROIDS] = clusterizeLockers(X, K) partitions the n-by-d data
%   matrix X into K clusters. IDX is an n-by-1 vector of cluster labels
%   in 1:K. CENTROIDS is K-by-d.
%
%   The algorithm enforces cluster size quotas so that each cluster gets
%   either floor(n/K) or ceil(n/K) points (difference at most 1).
%   It attempts to keep clusters compact by assigning points in order of
%   proximity to centroids at each iteration.
%
%   Optional name-value pairs:
%     'MaxIter'  - maximum iterations (default 100)
%     'Verbose'  - true/false (default false)
%     'Seed'     - random seed for reproducibility (default: rng untouched)
%
%   Example:
%     X = randn(300,2);
%     [idx,C] = clusterizeLockers(X, 5, 'MaxIter', 50, 'Seed',42);
%     gscatter(X(:,1), X(:,2), idx);
%     hold on; plot(C(:,1), C(:,2), 'kx','MarkerSize',12,'LineWidth',2);
%
%   Notes:
%     - k must be <= n (number of points).
%     - This is a heuristic algorithm (fast, robust) — it does not guarantee
%       global optimum but typically yields compact balanced clusters.
%

% Parse inputs
p = inputParser;
addRequired(p, 'X', @(x) isnumeric(x) && ndims(x)==2);
addRequired(p, 'k', @(x) isscalar(x) && x==floor(x) && x>0);
addParameter(p, 'MaxIter', 100, @(x) isscalar(x) && x>0);
addParameter(p, 'Verbose', false, @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'Seed', [], @(x) isempty(x) || (isscalar(x) && x==floor(x)));
parse(p, X, k, varargin{:});
maxIter = p.Results.MaxIter;
verbose = p.Results.Verbose;
seed = p.Results.Seed;

if ~isempty(seed)
    rng(seed);
end

[n, d] = size(X);
if k > n
    error('k (%d) must be <= number of points n (%d).', k, n);
end

% Compute quotas: each cluster gets q or q+1 points
q = floor(n / k);
r = mod(n, k); % number of clusters that should get q+1
quotas = repmat(q, k, 1);
if r > 0
    % to spread the extra ones, assign +1 to the first r clusters (we'll
    % reorder later by centroid position if desired)
    quotas(1:r) = quotas(1:r) + 1;
end
if sum(quotas) ~= n
    error('Internal quota error: quotas must sum to n.');
end

% --- kmeans++ initialization for centroids
centroids = zeros(k, d);
% choose first centroid uniformly at random from X
firstIdx = randi(n);
centroids(1, :) = X(firstIdx, :);

% distances to nearest chosen centroid (squared)
D2 = sum((X - centroids(1,:)).^2, 2);

for j = 2:k
    % choose next centroid with probability proportional to D2
    prob = D2 / sum(D2);
    cumprob = cumsum(prob);
    rdraw = rand();
    next = find(cumprob >= rdraw, 1, 'first');
    centroids(j, :) = X(next, :);
    % update D2
    newD2 = sum((X - centroids(j,:)).^2, 2);
    D2 = min(D2, newD2);
end

% main loop
idx = zeros(n,1);
prev_idx = idx;
for it = 1:maxIter
    % compute squared distances n x k
    % using bsxfun-like vectorization
    % D(i,j) = ||X(i,:) - centroids(j,:)||^2
    % Efficient compute:
    Cnorm = sum(centroids.^2, 2)';    % 1 x k
    Xnorm = sum(X.^2, 2);             % n x 1
    D = bsxfun(@plus, Cnorm, bsxfun(@minus, Xnorm, 2 * (X * centroids'))); % n x k
    % numeric safety
    D = max(D, 0);

    % For each point, get permutation of centroids sorted by distance
    [~, sortedIdxs] = sort(D, 2, 'ascend'); % n x k, each row: ordered centroids

    % order points by their minimum distance (closest-first)
    [minD, ~] = min(D, [], 2);
    [~, pointOrder] = sort(minD, 'ascend'); % indices 1..n

    % assignment respecting quotas
    remaining = quotas(:); % k x 1
    new_idx = zeros(n,1);

    for p_i = 1:n
        pt = pointOrder(p_i);
        orderList = sortedIdxs(pt, :);
        assigned = false;
        for cand = 1:k
            cluster = orderList(cand);
            if remaining(cluster) > 0
                new_idx(pt) = cluster;
                remaining(cluster) = remaining(cluster) - 1;
                assigned = true;
                break;
            end
        end
        if ~assigned
            % should not happen because sum remaining equals #unassigned,
            % but as a safeguard assign to any cluster with remaining>0
            cluster = find(remaining>0, 1, 'first');
            if isempty(cluster)
                % fallback: assign to nearest centroid
                new_idx(pt) = sortedIdxs(pt,1);
            else
                new_idx(pt) = cluster;
                remaining(cluster) = remaining(cluster) - 1;
            end
        end
    end

    % recompute centroids
    centroids_old = centroids;
    for j = 1:k
        members = (new_idx == j);
        if any(members)
            centroids(j, :) = mean(X(members, :), 1);
        else
            % This should not happen if quotas were used correctly, but in case:
            % reinitialize empty centroid to a random point
            centroids(j, :) = X(randi(n), :);
        end
    end

    % check for convergence
    if all(new_idx == idx)
        if verbose
            fprintf('balanced_kmeans: converged at iteration %d\n', it);
        end
        break;
    end

    idx = new_idx;

    if verbose
        % compute objective = sum squared distance to assigned centroid
        obj = sum( sum((X - centroids(idx,:)).^2, 2) );
        fprintf('it %d, objective %.6g\n', it, obj);
    end
end

% final outputs
idx = idx(:);
% ensure centroids correspond to final labels (they do)
end