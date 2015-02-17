function  df  = deriv( f, q )
%DERIV Computes the numerical derivative of a function f (passed by handle)
%
%   This function computes the numerical derivative of a function f at
%   position q. The function f is passed by handle as the following example
%   explains:
%
%   f  = @(x) sin(x);
%   df = deriv(f,0); 

[m , n] = size(q);
if n ~= 1
   error('In computing df(x), x should be a column vector')
end

[h , k] = size(f(q));
if k ~= 1 && m ~= 1
   error('In computing df(x), f should be a column vector or (if not) x should be scalar')
end

if m == 1
   q1 = q;
   h  = sqrt(eps)*max([abs(q1), 1]);
   q2 = q1 + h;
   
   f1 = f( q1 );
   f2 = f( q2 );
   
   df = (f2 - f1) ./ (q2 - q1);
   
else
   df = zeros(h,m);   
   for i = 1 : m
      q1 =  q;
      q2 =  q;
      h  = sqrt(eps)*max([abs(q1(i,1)), 1]);
      q2(i,1) = q1(i,1) + h;
      
      f1 = f( q1 );
      f2 = f( q2 );
      
      df(:,i) = (f2 - f1) ./ (q2(i,1) - q1(i,1));
   end
end


