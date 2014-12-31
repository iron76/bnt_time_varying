function  X = rotzderiv( theta )

% rotzderiv  derivative of spatial coordinate transform (Z-axis rotation).
% rotzderiv(theta)  calculates the derivative of function rotz with respect
% to theta

c = cos(theta);
s = sin(theta);

X = [ -s  c  0  0  0  0 ;
      -c -s  0  0  0  0 ;
       0  0  0  0  0  0 ;
       0  0  0 -s  c  0 ;
       0  0  0 -c -s  0 ;
       0  0  0  0  0  0
    ];