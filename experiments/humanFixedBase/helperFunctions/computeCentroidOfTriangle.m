function [ pCentroid ] = computeCentroidOfTriangle( p1, p2, p3 )
%COMPUTECENTROIDOFTRIANGLE Computes position of centroid of a trianlge by
%averaging
%   Assumes all points are expressed in some common reference

    pCentroid = (1/3)* ( p1+ p2+ p3);

end

