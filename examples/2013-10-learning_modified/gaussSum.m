clc
close
clear

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
bnet = buildGaussSum([1 1]);

engine{1} = pearl_inf_engine(bnet);

x1_m = 1;    x3_m = 3;    x2_m = 1;
evidence = cell(1,N);
evidence{x1} = x1_m;
evidence{x3} = x3_m;
[engine, ll] = enter_evidence(engine{1}, evidence);

s1 = 1;  m1  = 0;
s2 = 1;  m2  = 0;
s3 = 2;  m3  = 0;

marg = marginal_nodes(engine, x2)
s2/(s3+s2)*(x3_m - x1_m)
s3*s2/(s3+s2)


