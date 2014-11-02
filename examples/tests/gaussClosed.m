clc
close
clear

N = 4;
x1 = 1;  x2 = 3;  x3 = 4;  x4 = 2;
labels = {'x1', 'x2', 'x3', 'x4'};
dag = zeros(N,N);
dag(x1,x3) = 1;   dag(x2,x3) = 1;                     % x3 = x1 + x2
dag(x1,x4) = 1;   dag(x2,x4) = 1;   dag(x3,x4) = 1;   % x4 = x1 - x2 + x3


ns = ones(N,1); % no vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

x3_m = 1;   x4_m = 2;
sx1  = 1;   mx1  = 0;
sx2  = 4;   mx2  = 0;
sx3  = .02; mx3  = 0;
sx4  = .02; mx4  = 0;

Wx3     = [1 1];
Wx4     = [1 -1 1];
[~,ind] = sort([x1 x2]);
Wx3     = [Wx3(ind(1)) Wx3(ind(2))];
[~,ind] = sort([x1 x2 x3]);
Wx4     = [Wx4(ind(1)) Wx4(ind(2)) Wx4(ind(3))];

bnet.CPD{x1} = gaussian_CPD(bnet, x1, 'mean', mx1, 'cov', sx1);
bnet.CPD{x2} = gaussian_CPD(bnet, x2, 'mean', mx2, 'cov', sx2);
bnet.CPD{x3} = gaussian_CPD(bnet, x3, 'mean', mx3, 'cov', sx3, 'weights', Wx3);
bnet.CPD{x4} = gaussian_CPD(bnet, x4, 'mean', mx4, 'cov', sx4, 'weights', Wx4);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);

evidence{x3} = x3_m;
evidence{x4} = x4_m;

[engine, ll] = enter_evidence(engine, evidence);

marg_x1 = marginal_nodes(engine, x1);
marg_x1.mu
marg_x2 = marginal_nodes(engine, x2);
marg_x2.mu
draw_graph(bnet.dag, labels);
