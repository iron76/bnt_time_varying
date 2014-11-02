clc
close
clear

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
dag = zeros(N,N);
dag(1,3) = 1;
dag(2,3) = 1;

ns = [2 2 2]; % no vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

s1  = .1;  m1  = 0;
s2  = 10;  m2  = 0;
s12 =  2; m12  = 0;

m = [m1 m2 m12];
s = [s1 s2 s12];

for i=1:N
    if i == 3
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', [m(i); m(i)], 'cov', s(i).*eye(2), 'weights', [1 0 1 0; 0 1 0 1]);
    else
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', [m(i); m(i)], 'cov', s(i).*eye(2));
    end
end

engine{1} = gaussian_inf_engine(bnet);

evidence = cell(1,N);
x1_m = 1;    x3_m = 3;
evidence{x1} = [x1_m; x1_m];
evidence{x3} = [x3_m; x3_m];
[engine, ll] = enter_evidence(engine{1}, evidence);

marg = marginal_nodes(engine, x2)
s2/(s12+s2)*(x3_m - x1_m)
s12*s2/(s12+s2)