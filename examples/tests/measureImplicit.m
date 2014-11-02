clc
%%%%%%%%%%%%%%%%%%
% FIRST  EXAMPLE %
%%%%%%%%%%%%%%%%%%

clear all

N = 2;
h = 1;  k = 2;
dag = zeros(N,N);
dag(1,2) = 1;

ns = [3 4];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

mh = [ 0 0 0]';
Sk = [ .1 0 0 0; 0 .1 0 0 ; 0 0 .1 0; 0 0 0 .1];
dumpFactor = 1e4;
Sh  = dumpFactor.*[ 1 0 0; 0 1 0; 0 0 1];

Ch = [ 1 0  0; 0 1 0; 0 0 1; 0 0 0];
Cy = [ 1 1; 2 1; 3 -1; 1 1];

bnet.CPD{h} = gaussian_CPD(bnet, h, 'mean', mh, 'cov', Sh);
bnet.CPD{k} = gaussian_CPD(bnet, k, 'mean', [0 0 0 0]', 'cov', Sk, 'weights', Ch);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
evidence{k} = Cy * [1 1]';
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, h);
Ehk = marg.mu;
Shk = marg.Sigma;

% It results: 
% Ehk
% (Sh^(-1) + Ch'*Sk^(-1)*Ch)^(-1)*Ch'*Sk^(-1)*Cy*[1 1]'
% Shk
% (Ch'*Sk^(-1)*Ch)^(-1)


N = 3;
y = 1;  ym = 2;  h = 3;
dag = zeros(N,N);
dag(1,2) = 1;    dag(1,3) = 1;

ns = [2 2 3];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% Cx x + Cy y = w
%
% x ~ N(mx, Sx)
% w ~ N(0 , R)

my   = [ 0 0]';
Sy   = [ 2 1; 1 2];
Sym  = [.1 0; 0 .1];

bnet.CPD{y}  = gaussian_CPD(bnet, y,  'mean', my, 'cov', Sy);
bnet.CPD{ym} = gaussian_CPD(bnet, ym, 'mean', [0 0]', 'cov', Sym, 'weights', eye(2));
bnet.CPD{h}  = gaussian_CPD(bnet, h,  'mean', [0 0 0]', 'cov', Shk, 'weights', -Shk*Ch'*Sk^(-1)*Cy);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
evidence{ym} = [1 1]';
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [y h]);
Ehm = marg.mu;
Shm = marg.Sigma;

% It results: 
% Ehm
% (Ch'*Sk^(-1)*Ch)^(-1)*Ch'*Sk^(-1)*Cy*[1 1]'
Shm
(Ch'*(Sk + Cy * Sym*(Sy+Sym)^(-1)*Sy*Cy')^(-1) * Ch)^(-1)

%%%%%%%%%%%%%%%%%%%
% SECOND  EXAMPLE %
%%%%%%%%%%%%%%%%%%%
clear all

N = 3;
y = 1;  ym = 2;  k= 3;
dag = zeros(N,N);
dag(y,ym) = 1;    dag(y,k) = 1;

ns = [2 2 4];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% Cy x + Ch h = w
% Ch h = k
% ym = y + v
% 
%
% x ~ N(mx, Sx)
% w ~ N(0 , Sk)
% v ~ N(0, Sym)

my   = [ 0 0]';
Sy   = [ 2 1; 1 2];
Sym  = [.1 0; 0 .1];

Sk = [ .1 0 0 0; 0 .1 0 0 ; 0 0 .1 0; 0 0 0 .1];

Ch = [ 1 0  0; 0 1 0; 0 0 1; 0 0 0];
Cy = [ 1 1; 2 1; 3 -1; 1 1];


bnet.CPD{y}  = gaussian_CPD(bnet, y,  'mean', my, 'cov', Sy);
bnet.CPD{ym} = gaussian_CPD(bnet, ym, 'mean', [0 0]', 'cov', Sym, 'weights', eye(2));
bnet.CPD{k}  = gaussian_CPD(bnet, k,  'mean', [0 0 0 0]', 'cov', Sk, 'weights', Cy);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
evidence{ym} = [1 1]';
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [y k]);
Ekm = marg.mu;
Skm = marg.Sigma;

% It results: 
Skm
% (Sk + Cy * Sym*(Sy+Sym)^(-1)*Sy*Cy')

% If necessary we can also retrieve the estimation of h
(Ch' * Skm(3:6, 3:6)^(-1) * Ch)^(-1)
% and the the covariance of y-h (to be done)


