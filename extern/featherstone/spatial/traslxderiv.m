function  X = traslxderiv( theta )

% rotxderiv  derivative of spatial coordinate transform (X-axis rotation).
% rotxderiv(theta)  calculates the derivative of function rotx with respect
% to theta


X = [  0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0     0    0  0  0 ;
       0     0     1    0  0  0 ;
       0    -1     0    0  0  0
    ];
