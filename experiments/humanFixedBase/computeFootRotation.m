function [ R_G_0 , P_G_0] = computeFootRotation( P_G_lheeT,P_G_rheeT, P_G_ltoeT, P_G_rtoeT)
%COMPUTEFOOTROTATION Computes the elements of the roation matrix of the
%foot frame, frame in link 0
%   Detailed explanation goes here
    numPoints = 10;
    P_G_lhee = mean(P_G_lheeT(1:numPoints,:));
    P_G_rhee = mean(P_G_rheeT(1:numPoints,:));
    P_G_ltoe = mean(P_G_ltoeT(1:numPoints,:));
    P_G_rtoe = mean(P_G_rtoeT(1:numPoints,:));
    
    
    P_G_mhee = computeCentroidOfPoints(P_G_lhee,P_G_rhee);
    P_G_mtoe = computeCentroidOfPoints(P_G_ltoe,P_G_rtoe);    
    P_G_0 = computeCentroidOfPoints(P_G_mhee,P_G_mtoe);
    
    v_G_0 = computeVectorFromPoints(P_G_0,P_G_mtoe);
    x_G_0 = v_G_0 ./ norm(v_G_0);
    z_G_0 = repmat([0,0,1],size(v_G_0,1),1);
    y_G_0 = cross(x_G_0,z_G_0);
    
    R_G_0 = [x_G_0', y_G_0', z_G_0']; 

end

