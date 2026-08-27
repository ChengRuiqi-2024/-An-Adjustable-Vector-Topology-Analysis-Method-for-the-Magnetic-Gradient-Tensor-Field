function result = classify_critical_points_rectangular_ring_with_focus_candidates(X,Y,U,V,rowRadius,colRadius,opts)
% Discrete critical-point classification on a user-defined rectangular ring.
%
% INPUT
%   X,Y,U,V       equally sized grid matrices (as produced by meshgrid)
%   rowRadius     centre-to-boundary distance in grid intervals vertically
%   colRadius     centre-to-boundary distance in grid intervals horizontally
%
% Thus rowRadius = 1, colRadius = 1 reproduces the original 3-by-3
% intervals and contains the eight surrounding points.
% result.type:  1 source-like; -1 sink-like; 2 saddle-like;
%               3 non-radial +1-index candidate (including focus-like points);
%               0 unclassified
% result.index: discrete Poincare index on the rectangular ring
% result.score: mean signed radial-direction cosine on the ring
% result.nRingPoints: number of boundary samples actually used
% result.validRing: true where a complete, nonzero and unambiguous ring exists
%
% This method does not interpolate and does not calculate a Jacobian.
% Direction-only classification identifies local topology candidates;
arguments
    X double
    Y double
    U double
    V double
    rowRadius (1,1) double {mustBeInteger,mustBePositive}
    colRadius (1,1) double {mustBeInteger,mustBePositive}
    opts.radialCosThreshold (1,1) double = 0.30
    opts.minConsistentFraction (1,1) double = 0.875
    opts.indexTolerance (1,1) double = 0.55
    opts.maxDirectionStep (1,1) double = 0.98*pi
    opts.centerMagnitudeThreshold (1,1) double = inf
end

[nr,nc] = size(U);
assert(isequal(size(X),[nr nc],size(Y),[nr nc],size(V),[nr nc]), ...
    'X, Y, U and V must have identical sizes.');

% All grid points on the rectangular boundary, expressed relative to its centre.
[dcol,drow] = meshgrid(-colRadius:colRadius, -rowRadius:rowRadius);
isBoundary = abs(drow)==rowRadius | abs(dcol)==colRadius;
drow = drow(isBoundary);
dcol = dcol(isBoundary);

result.type = zeros(nr,nc,'int8');
result.index = nan(nr,nc);
result.score = nan(nr,nc);
result.nRingPoints = zeros(nr,nc,'uint16');
result.validRing = false(nr,nc);

for i = 1+rowRadius:nr-rowRadius
    for j = 1+colRadius:nc-colRadius
        ii = i + drow;
        jj = j + dcol;
        lin = ii + (jj-1)*nr;

        xr = X(lin);  yr = Y(lin);
        ur = U(lin);  vr = V(lin);
        speed = hypot(ur,vr);

        % A Poincare index is undefined if a boundary vector has no direction.
        if any(~isfinite(xr)) || any(~isfinite(yr)) || ...
           any(~isfinite(speed)) || any(speed == 0)
            continue
        end

        % Sort actual physical positions counter-clockwise.  This makes the
        % result independent of whether rows of Y increase upward or downward.
        alpha = atan2(yr-Y(i,j), xr-X(i,j));
        [~,order] = sort(alpha,'ascend');
        xr = xr(order);  yr = yr(order);
        ur = ur(order);  vr = vr(order);  speed = speed(order);

        phi = atan2(vr,ur);
        dphi = atan2(sin(diff([phi; phi(1)])), cos(diff([phi; phi(1)])));

        % If an adjacent sampled direction turns nearly pi, its true rotation
        % is ambiguous without interpolation; leave this centre unclassified.
        if any(abs(dphi) > opts.maxDirectionStep)
            continue
        end

        ind = sum(dphi)/(2*pi);
        rx = xr-X(i,j);
        ry = yr-Y(i,j);
        radialCos = (ur.*rx + vr.*ry) ./ (speed .* hypot(rx,ry));

        result.index(i,j) = ind;
        result.score(i,j) = mean(radialCos);
        result.nRingPoints(i,j) = numel(ur);
        result.validRing(i,j) = true;

        % Optional magnitude gate; leave Inf to use directions only.
        if hypot(U(i,j),V(i,j)) > opts.centerMagnitudeThreshold
            continue
        end

        fracOut = mean(radialCos >  opts.radialCosThreshold);
        fracIn  = mean(radialCos < -opts.radialCosThreshold);

        if abs(ind - 1) <= opts.indexTolerance && ...
                fracOut >= opts.minConsistentFraction
            result.type(i,j) = 1;
        elseif abs(ind - 1) <= opts.indexTolerance && ...
                fracIn >= opts.minConsistentFraction
            result.type(i,j) = -1;
        elseif abs(ind - 1) <= opts.indexTolerance
            % Index +1, but neither predominantly inward nor outward:
            % non-radial candidate, including noise-induced focus-like points.
            result.type(i,j) = 3;
        elseif abs(ind + 1) <= opts.indexTolerance
            result.type(i,j) = 2;
        end
    end
end

end
