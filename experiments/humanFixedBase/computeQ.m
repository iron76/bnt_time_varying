function [Q] = computeQ(map, q,dq,Y,y,b_y,theta1,theta2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q,dq,Y, y,b_y, inv(theta1));
    yHat = y - b_y;
    
    n_y = map.IDsens.m;
    n_D = length(map.d);


    Q1 = -0.5*(n_y+n_D)*log(2*pi) + 0.5*log(det(inv(map.sigmaBarD))) + 0.5*log(det(inv(theta2))); 
    Q2 = -0.5*(mu_dgiveny-map.muBarD)'* inv(map.sigmaBarD)* (mu_dgiveny-map.muBarD) - 0.5*trace(inv(map.sigmaBarD) * Sigma_dgiveny);
    Q3 = -0.5*((Y*mu_dgiveny)-yHat)'*inv(theta2)*((Y*mu_dgiveny)-yHat) - 0.5*trace(inv(theta2)*Y*Sigma_dgiveny*Y');
    Q =  (Q1 + Q2 + Q3);
   
    
% %         [mu_dgiveny, Sigma_dgiveny, map] = computeMAP(map,q,dq,Y, y,b_y, inv(theta1));
% %     yHat = y - b_y;
% %     
% %     n_y = map.IDsens.m;
% %     n_D = length(map.d);
% % 
% % 
% %     Q1 = -0.5*(n_y+n_D)*log(2*pi) + 0.5*log(det(inv(map.sigmaBarD))) + 0.5*log(det(inv(theta2))); 
% %     Q2 = -0.5*map.muBarD'*(inv(map.sigmaBarD))*map.muBarD - 0.5*yHat'*(inv(theta2))*yHat + map.muBarD'*(inv(map.sigmaBarD))*mu_dgiveny + mu_dgiveny'*Y'*(inv(theta2))*yHat;
% %     Q3 = -0.5*trace((inv(map.sigmaBarD)+(Y'*inv(theta2)*Y))*(Sigma_dgiveny)) - 0.5*map.muBarD'*(inv(map.sigmaBarD)+(Y'*inv(theta2)*Y))*map.muBarD; 
% %     Q =  (Q1 + Q2 + Q3);
% %     
end

