% % % % function [l_y] = computeLogL(map, q,dq,Y,y,b_y,theta)
% % % % %UNTITLED2 Summary of this function goes here
% % % % %   Detailed explanation goes here
% % % % 
% % % %     [~, ~, map] = computeMAP(map,q,dq,Y, y,b_y, inv(theta));
% % % %     yHat = y - b_y;
% % % % 
% % % %     n_y = map.IDsens.m;
% % % %     SigmaD = map.sigmaBarD;
% % % %     c = yHat-(Y*map.muBarD);
% % % %     Sigma_y = theta+(Y*SigmaD*Y');
% % % %     detSigma = det(Sigma_y);
% % % % 
% % % %     l_y = -0.5*n_y*log(2*pi) -0.5*log(detSigma) -0.5*(c'*inv(Sigma_y)*c);
% % % %   
% % % % end


function [l1] = computeLogL(map, q,dq,Y,y,b_y,theta1,theta2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    [~, ~, map] = computeMAP(map,q,dq,Y, y,b_y, inv(theta1));
    yHat = y - b_y;

    n_y = map.IDsens.m;
    n_D = length(map.d);
    
    l1 = -0.5*(n_D+n_y)*log(2*pi) +0.5*log(det(inv(map.sigmaBarD)))+0.5*log(det(inv(theta2)))...
         -0.5*(map.d-map.muBarD)'*inv(map.sigmaBarD)*(map.d-map.muBarD) ...
         -0.5*(yHat-(Y*map.muBarD))'*inv(theta2)*(yHat-(Y*map.muBarD));
    
% %     p_d =  ((2*pi)^(-0.5*n_D ) * ((det(theta1))^(-0.5)) * exp( -0.5*(map.d-map.muBarD)'*inv(map.sigmaBarD)*(map.d-map.muBarD)) );
% %     p_ygivend =  ((2*pi)^(-0.5*n_y ) * ((det(theta2))^(-0.5)) * exp( -0.5*(yHat-(Y*map.muBarD))'*inv(theta2)*(yHat-(Y*map.muBarD))) );
% %     
% %     
% %    l1 = log(p_d*p_ygivend);
   
end


% function [l1,l2] = computeLogL(map, q,dq,Y,y,b_y,theta)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% 
%     [~, ~, map] = computeMAP(map,q,dq,Y, y,b_y, inv(theta));
%     yHat = y - b_y;
% 
%     n_y = map.IDsens.m;
%     
%     l1 =(-0.5*(n_y*log(2*pi)+log(det(theta))))- (0.5*(yHat-(Y*map.muBarD))'*inv(theta)*(yHat-(Y*map.muBarD)) );
%     l2 = log( ((2*pi)^(-0.5*n_y ) * ((det(theta))^(-0.5)) * exp( -0.5*(yHat-(Y*map.muBarD))'*inv(theta)*(yHat-(Y*map.muBarD))) ));
% end
