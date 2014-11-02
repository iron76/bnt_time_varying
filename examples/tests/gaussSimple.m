clc
close
clear

N = 2;
x = 1;  y = 2;
dag = zeros(N,N);
dag(x,y) = 1;

ns = [3 2]; % no vector-valued
S = [1 2 3; 1 1 0];

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

mx = [1 0 0]';
my = [2 1]';
Sx = [3 1 1; 0 2 1; 0 0 1];
Sx = Sx + Sx';
Sy = [2 1; 1 2];
Sx = Sx + Sx';
ym = [1 1]';

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sx);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', my, 'cov', Sy, 'weights', S);   % y = Sx + w

engine = jtree_inf_engine(bnet);

evidence = cell(1,N);
evidence{y} = ym;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
marg.mu
mx + (Sx^(-1) + S'*Sy^(-1)*S)^(-1)*S'*Sy^(-1)*((ym-my)-S*mx)

marg.Sigma
(Sx^(-1) + S'*Sy^(-1)*S)^(-1)


