function [ a_filt ] = timeFilterSignal( polynOrder,window, a_nofilt)

% TIMEFILTERSIGNAL computes the filtering of a signal that has to be in the
% form (lengthSignal,3).  As it uses SGolay Smoother, it requires also the 
% polynomial order and the size of the moving window.

%a_filt = a_nofilt;
a_filt = zeros(length(a_nofilt),3);
    for i = 1:3
         a_filt(:,i) = SgolaySmoother(polynOrder,window,a_nofilt(:,i));
    end
end
