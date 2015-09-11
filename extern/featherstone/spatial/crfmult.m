function [ vcross ] = crfmult( v, f )
%crfmult computes crf(v)*f


if length(v) == 6

%   vcross = -[  0    -v(3)  v(2)   0     0     0    ;
% 	      v(3)  0    -v(1)   0     0     0    ;
% 	     -v(2)  v(1)  0      0     0     0    ;
% 	      0    -v(6)  v(5)   0    -v(3)  v(2) ;
% 	      v(6)  0    -v(4)   v(3)  0    -v(1) ;
% 	     -v(5)  v(4)  0     -v(2)  v(1)  0 ]'*f;
     
%   vcross =  [  0    -v(3)  v(2)   0    -v(6)  v(5)    ;
% 	              v(3)  0    -v(1)   v(6)  0    -v(4)    ;
% 	             -v(2)  v(1)  0     -v(5)  v(4)  0    ;
% 	              0     0     0      0    -v(3)  v(2) ;
% 	              0     0     0      v(3)  0    -v(1) ;
% 	              0     0     0     -v(2)  v(1)  0 ]*f;
     
  vcross =  [  -v(3)*f(2)  + v(2)*f(3) - v(6)*f(5)  + v(5)*f(6)    ;
	              v(3)*f(1) - v(1)*f(3) + v(6)*f(4)  - v(4)*f(6)    ;
	             -v(2)*f(1) + v(1)*f(2) - v(5)*f(4)  + v(4)*f(5)    ;
	              -v(3)*f(5)+ v(2)*f(6) ;
	              v(3)*f(4) - v(1)*f(6) ;
	              -v(2)*f(4) + v(1)*f(5)];     
     % vcross = [cross(v(1:3,1), f(1:3,1)) + cross(v(4:end,1), f(4:end,1)); cross(v(1:3,1), f(4:end,1))];

else

  vcross = -[  0     0     0    ;
	      v(3)  0    -v(1) ;
	     -v(2)  v(1)  0 ]'*f;
end

end

