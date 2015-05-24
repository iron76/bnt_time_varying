function [ pCentroid ] = computeCentroidOfPoints( p1,p2 )
%COMPUTECENTROIDOFPOINTS Computes the location of the centroid of 2 points
%   Points expressed in some common reference frame

pCentroid = 0.5*(p1+p2);

end

