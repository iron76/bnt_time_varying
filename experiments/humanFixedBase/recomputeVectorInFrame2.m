function [P_2_X] = recomputeVectorInFrame2(q1,R_1_G,P_G_X)
        
    P_2_X = zeros(size(P_G_X));
    %R_1_2 = cell(size(P_G_X,1));
    P_1_X = (R_1_G*P_G_X')';
    for i = 1:length(q1)
   %     R_1_2{i} = euler2dcm(q1(i),0,0);
        R_1_2 = euler2dcm([q1(i),0,0]);
        P_2_X(i,:) = R_1_2'*P_1_X(i,:)';
    end
end