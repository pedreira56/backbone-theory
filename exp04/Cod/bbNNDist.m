function [R, J] = bbNNDist(Xq, Xd, P)
%BBNNDIST  Distance from each receiver to its nearest donor in backbone space.
%
%   R = BBNNDIST(XQ,XD,P) returns, for each row of the nq-by-b matrix XQ of
%   receiver backbones, the distance to the nearest row of the nd-by-b matrix
%   XD of donor backbones.  [R,J] = ... also returns the donor indices.
%
%   The metric follows P.geometry:
%     'cube'   Euclidean, evaluated with a kd-tree.
%     'torus'  flat torus on [0,1)^d with wrap-around in every coordinate,
%              evaluated exactly in chunks.  The torus has no boundary, so
%              Cond_G holds with a uniform constant at every point of the
%              support; it is the clean control for E1, against which the
%              cube shows the boundary-inflated constant at the same rate.
%
%   See also BBSAMPLE, E1_NNRATE.

if size(Xq,2) ~= size(Xd,2)
    error('bbNNDist:dimMismatch','XQ and XD must have the same number of columns.');
end

switch P.geometry
    case 'cube'
        ns = createns(Xd,'NSMethod','kdtree');
        [J,R] = knnsearch(ns,Xq,'K',1);

    case 'torus'
        nq = size(Xq,1);
        R  = zeros(nq,1);
        J  = zeros(nq,1);
        % Chunk the receivers so the nq-by-nd-by-d difference array stays
        % within a few hundred MB regardless of the sample size.
        bytesPerChunk = 200e6;
        chunk = max(1, floor(bytesPerChunk / (8 * size(Xd,1) * size(Xd,2))));
        for i0 = 1:chunk:nq
            i1 = min(nq, i0+chunk-1);
            A  = reshape(Xq(i0:i1,:), [], 1, size(Xq,2));
            B  = reshape(Xd,          1, [], size(Xd,2));
            D  = abs(A - B);
            D  = min(D, 1 - D);                 % wrap-around per coordinate
            D  = sqrt(sum(D.^2, 3));            % (i1-i0+1) x nd
            [rv,jv] = min(D, [], 2);
            R(i0:i1) = rv;
            J(i0:i1) = jv;
        end

    otherwise
        error('bbNNDist:badGeometry','Unknown geometry "%s".',P.geometry);
end
end
