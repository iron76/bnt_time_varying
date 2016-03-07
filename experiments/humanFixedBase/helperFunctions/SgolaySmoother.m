function [ a_smoot ] = SgolaySmoother( polynOrder,window,a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%diffCoeff is a matrix of (polynomialOrder-1) columns where:
%- ( ,1) --> coefficient for S-Golay as smoother;
%- ( ,2) --> coefficient for S-Golay as 1st differentiator;
%- ( ,3) --> coefficient for S-Golay as 2nd differentiator;
%  .   
%  .
%  .
%- ( ,polynomialOrder-1) --> coefficient for S-Golay as (polynomialOrder) differentiator;




[~, diffCoeff] = sgolay_wrapper(polynOrder, window);
halfWindow  = ((window+1)/2) -1;
l = length(a);
a_smoot = zeros(l, 1);
 
 for n = (window+1)/2:l-(window+1)/2,        
     a_smoot(n) = dot(diffCoeff(:,1),a(n - halfWindow:n + halfWindow));             
 end

end









