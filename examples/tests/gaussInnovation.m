clc
close
clear

% We here prove that also in cases which the distribution of y is 
% implicitly defined as y = C x + w with x ~ N(mx, Sx) and
% w ~ N(0 , R) we can compute E[x|y] with a graphical model regardless
% of the dimensions of x and y as long as Sx and R are invertible.

%%%%%%%%%%%%%%%%
% FIRST EXAMPLE%
%%%%%%%%%%%%%%%%
N = 2;
x    = 1;     y = 2;    %variables ID
s(x) = 3;  s(y) = 4;    %variables sizes
dag = zeros(N,N);       %connection matrix
dag(x,y) = 1;


dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, s, 'discrete', dnodes);

ym = [ 1 1 1 1]';
my = [ 0 0 0 0]';
mx = [ 1 0 -1 ]';
Sx = [ 1 0 -1/2; 0 1 0; -1/2 0 1];
R  = [ 1 0 0 0; 0 2 0 0; 0 0 1 0 ; 0 0 1 1];

C  = [ 1 2 3; 4 5 6; 0 0 1; 1 1 1];

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sx);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', my, 'cov', R, 'weights', C);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
evidence{y} = ym;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
Exy = marg.mu;
Sxy = marg.Sigma;

% It results: 
% Exy = mx + Sx*C'*(R+C*Sx*C')^(-1)*(ym - C*mx)
% Sxy = Sx - Sx*C'*(R+C*Sx*C')^(-1)*C*Sx 

%%
%%%%%%%%%%%%%%%%%%
% SECOND EXAMPLE % 
%%%%%%%%%%%%%%%%%%

% We here solve the same problem as before, but in this case we separate 
% estimation in two steps. The model y = C x + w with x ~ N(mx, Sx) 
% with w ~ N(0 , R) is divided in two steps, expliting the fact that the
% matrix R has a block diagonal structure. In particular, we 
% R = diag(R1, R2), C = [C1; C2], y = [y1; y2] and therefore the estimation
% is divided in two equations y1 = C1 x + w1 and  y2 = C2 x + w2. The
% solution is obtained by first solving x | y1, which is then used as the
% prior for the second estimation based on y2 = C2 x + w2. 


%% First step

N = 2;
x    = 1;     y1 = 2;    %variables ID
s(x) = 3;  s(y1) = 2;    %variables sizes
dag = zeros(N,N);       %connection matrix
dag(x,y1) = 1;


dnodes = [];  % no discrete nodes
bnet1 = mk_bnet(dag, s, 'discrete', dnodes);


mx  = [ 1 0 -1 ]';
Sx  = [ 1 0 -1/2; 0 1 0; -1/2 0 1];
R1  = R(1:2,1:2); 
C1  =  C(1:2, :);
my1 = my(1:2, :); 
ym1 = ym(1:2, :); 

bnet1.CPD{x}  = gaussian_CPD(bnet1,  x, 'mean', mx,  'cov', Sx);
bnet1.CPD{y1} = gaussian_CPD(bnet1, y1, 'mean', my1, 'cov', R1, 'weights', C1);

engine = gaussian_inf_engine(bnet1);

evidence = cell(1,N);
evidence{y1} = ym1;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
Exy1 = marg.mu;
Sxy1 = marg.Sigma;

%% Second step

N = 2;
x    = 1;     y2 = 2;    %variables ID
s(x) = 3;  s(y2) = 2;    %variables sizes
dag = zeros(N,N);       %connection matrix
dag(x,y2) = 1;


dnodes = [];  % no discrete nodes
bnet2 = mk_bnet(dag, s, 'discrete', dnodes);


mx  = Exy1;
Sx  = Sxy1;
R2  = R(3:4,3:4); 
C2  =  C(3:4, :);
my2 = my(3:4, :); 
ym2 = ym(3:4, :); 

bnet2.CPD{x}  = gaussian_CPD(bnet2,  x, 'mean',  mx,  'cov', Sx);
bnet2.CPD{y2} = gaussian_CPD(bnet2, y2, 'mean', my2, 'cov', R2, 'weights', C2);

engine = gaussian_inf_engine(bnet2);

evidence = cell(1,N);
evidence{y2} = ym2;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
Exy2 = marg.mu;
Sxy2 = marg.Sigma;


