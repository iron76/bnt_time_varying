function res = testAnalytical_NumericalDerivative (f, x, df_analyt)

% Test to compare the numerical derivative of a function f and its 
% anlytical derivative (passed as an input).
% The function f is passed by handle as follows:
% f  = @(x) sin(x);

[m,n]  = size (x);
if m ~= 1 && n == 1
     error('In this test, x should be a scalar or a matrix')
end

rng(0);
res = 0;
   
df_num   = deriv(f, x);   

if (df_num - df_analyt) > 1e-8;
    disp(['[DERIVATIVES] Numerical and analytical derivative are different!']);
    res = 1;
else
    disp(['[DERIVATIVES] Numerical and analytical derivative are similar.']);
end

end

