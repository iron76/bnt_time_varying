function [ val ] = rowWiseNorm( v )
%ROWWISENORM Summary of this function goes here
%   Detailed explanation goes here

    val = sqrt(sum(v.^2,2));
end

