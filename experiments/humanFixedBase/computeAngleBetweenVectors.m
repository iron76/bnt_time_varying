function [ theta ] = computeAngleBetweenVectors( v1,v2 )
%COMPUTEANGLEBETWEENVECTORS compute the angle from v1 to v2
%   

  %  theta = acos(dot(v1',v2')'./...
  %      ((rowWiseNorm(v1) + rowWiseNorm(v2)).^0.5));

    
 %   theta = asin(rowWiseNorm(cross(v1,v2))./...
 %       (rowWiseNorm(v1).*rowWiseNorm(v2)));

    cv = cross(v1,v2);
    cvN = rowWiseNorm(cv);
   theta = atan2(cvN,dot(v1',v2')');

%    theta  = 2 * atan2(rowWiseNorm(v1.*rowWiseNorm(v2) - rowWiseNorm(v1).*v2) , rowWiseNorm(v1.*rowWiseNorm(v2) + rowWiseNorm(v1).*v2));
    
end

