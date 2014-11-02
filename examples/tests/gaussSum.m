clc
close
clear

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
dag = zeros(N,N);
dag(1,3) = 1;
dag(2,3) = 1;

ns = [1 1 1]; % no vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

x1_m = 1;    x3_m = 3;    x2_m = 1;
s1 = 1;  m1  = 0;
s2 = 1;  m2  = 0;
s3 = .02;  m3  = 0;

m = [m1 m2 m3];
s = [s1 s2 s3];

% for i=1:N
%     if i == 3
%         bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', m(i), 'cov', s(i), 'weights', [1 1]);   % x3 = x1 + x2
%     else
%         bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', m(i), 'cov', s(i));
%     end
% end

bnet.CPD{x1} = gaussian_CPD(bnet, x1, 'mean', m(x1), 'cov', s(x1));
bnet.CPD{x2} = gaussian_CPD(bnet, x2, 'mean', m(x2), 'cov', s(x2));
bnet.CPD{x3} = gaussian_CPD(bnet, x3, 'mean', m(x3), 'cov', s(x3), 'weights', [1 1]);

engine = var_elim_inf_engine(bnet);

evidence = cell(1,N);
evidence{x1} = x1_m;
evidence{x3} = x3_m;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x2)
s2/(s3+s2)*(x3_m - x1_m)
s3*s2/(s3+s2)

