function  df  = deriv( f, q )
%DERIV Computes the numerical derivative of a function f (passed by handle)
%
%   This function computes the numerical derivative of a function f at
%   position q. The function f is passed by handle as the following example
%   explains:
%
%   f  = @(x) sin(x);
%   df = deriv(f,0); 

q1 = q;
h  = sqrt(eps)*max([abs(q1), 1]);
q2 = q1 + h;

f1 = f( q1 );
f2 = f( q2 );

df = (f2 - f1) ./ (q2 - q1);

end

