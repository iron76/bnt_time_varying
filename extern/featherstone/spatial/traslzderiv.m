function  X = traslzderiv( theta )

% rotzderiv  derivative of spatial coordinate transform (Z-axis rotation).
% rotzderiv(theta)  calculates the derivative of function rotz with respect
% to theta


X = [  0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     1     0    0  0  0 ;
      -1     0     0    0  0  0 ;
       0     0     0    0  0  0
    ];