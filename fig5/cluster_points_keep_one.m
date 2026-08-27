function [representativeXY, clusterID, representativeIndex] = ...
    cluster_points_keep_one(xy, clusterRadius)
%  Perform spatial clustering on the two-dimensional coordinates, and retain one original point for each cluster.
%
% input：
%   xy            N×2 array, with one point coordinate [x, y] in each row
%   clusterRadius Clustering radius; the same as unit and xy
%
% output：
%   representativeXY     K×2, the coordinate of one representative point retained for each cluster
%   clusterID            N×1, the cluster number to which each original point belongs
%   representativeIndex  K×1, the row number of the representative point in the original xy coordinate system


     arguments
        xy (:,2) double
        clusterRadius (1,1) double {mustBePositive}
     end
        valid = all(isfinite(xy), 2);

    if ~any(valid)
        representativeXY = zeros(0,2);
        representativeIndex = zeros(0,1);
        return
    end

% minPts = 1: An isolated point itself constitutes a cluster
    clusterID = dbscan(xy, clusterRadius, 1);

    clusterLabels = unique(clusterID, 'stable');

    representativeIndex = zeros(numel(clusterLabels), 1);

    for k = 1:numel(clusterLabels)

        thisLabel = clusterLabels(k);

% Index of the original points in the current cluster
        idx = find(clusterID == thisLabel);

        pts = xy(idx, :);

% Select the "original point" that is closest to the average position of this cluster as the representative point.
        centre = mean(pts, 1);

        [~, pos] = min(sum((pts - centre).^2, 2));

        representativeIndex(k) = idx(pos);
    end

    representativeXY = xy(representativeIndex, :);

end