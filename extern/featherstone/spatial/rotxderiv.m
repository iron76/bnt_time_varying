function  X = rotxderiv( theta )

% rotxderiv  derivative of spatial coordinate transform (X-axis rotation).
% rotxderiv(theta)  calculates the derivative of function rotx with respect
% to theta

c = cos(theta);
s = sin(theta);

X = [ 0  0  0  0  0  0 ;
      0 -s  c  0  0  0 ;
      0 -c -s  0  0  0 ;
      0  0  0  0  0  0 ;
      0  0  0  0  s  c ;
      0  0  0  0 -c  s
    ];
