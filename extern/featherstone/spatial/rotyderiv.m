function  X = rotyderiv( theta )

% rotyderiv  derivative of spatial coordinate transform (Y-axis rotation).
% rotyderiv(theta)  calculates the derivative of function roty with respect
% to theta

c = cos(theta);
s = sin(theta);

X = [-s  0 -c  0  0  0 ;
      0  0  0  0  0  0 ;
      c  0 -s  0  0  0 ;
      0  0  0 -s  0 -c ;
      0  0  0  0  0  0 ;
      0  0  0  c  0 -s
    ];
