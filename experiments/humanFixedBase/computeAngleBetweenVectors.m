function [ theta ] = computeAngleBetweenVectors( v1,v2 )
%COMPUTEANGLEBETWEENVECTORS compute the angle from v1 to v2
%   

%    theta = acos(dot(v1',v2')'./...
%        (rowWiseNorm(v1).*rowWiseNorm(v2)));

    
    %theta = asin(rowWiseNorm(cross(v1,v2))./...
     %   (rowWiseNorm(v1).*rowWiseNorm(v2)));

    theta = atan2(rowWiseNorm(cross(v1,v2)),dot(v1',v2')');
end

