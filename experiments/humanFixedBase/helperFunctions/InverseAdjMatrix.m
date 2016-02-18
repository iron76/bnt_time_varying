function B = InverseAdjMatrix(A)

% InverseAdjMatrix provides the inverse transformation from 

% A_X_B = [    R      0          
%           -RS(r)    R ]               
% 
%    to
%
% B_X_A = [   R^T      0  
%           S(r)R^T   R^T ]


block1 = A(1:3,1:3); %R
block2 = A(1:3,4:6); %0
block3 = A(4:6,1:3); %-RS(r)
block4 = A(4:6,4:6); %R

B = [block1'  block2;
     block3'  block4'];
 
end

