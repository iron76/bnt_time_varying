function  X = traslyderiv( theta )

% rotyderiv  derivative of spatial coordinate transform (Y-axis rotation).
% rotyderiv(theta)  calculates the derivative of function roty with respect
% to theta


X = [  0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0    -1    0  0  0 ;
       0     0     0    0  0  0 ;
       1     0     0    0  0  0
    ];
