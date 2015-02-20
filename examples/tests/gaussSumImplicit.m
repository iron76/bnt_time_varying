clc
close
clear

% We here prove that also in cases which the distribution of y is 
% implicitly defined as Cy y + Cx x = w with x ~ N(mx, Sx) and
% w ~ N(0 , R) we can compute E[x|y] with a graphical model regardless
% of the dimensions of x and y as long as Sx and R are invertible.

%%%%%%%%%%%%%%%%
% FIRST EXAMPLE%
%%%%%%%%%%%%%%%%
N = 2;
x = 1;  y = 2;
dag = zeros(N,N);
dag(1,2) = 1;

ns = [3 3];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% Cx x + Cy y = w
%
% x ~ N(mx, Sx)
% w ~ N(0 , R)

mx = [ 1 0 -1 ]';
Sx = [ 1 0 -1/2; 0 1 0; -1/2 0 1];
R  = [ 1 0 0; 0 2 0; 0 0 1];

Cx = [ 1 2 3; 4 5 6; 0 0 1];
Cy = [ 1 2; 2 3; 1 1];

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sx);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', [0 0 0]', 'cov', R, 'weights', -Cx);

engine{1} = gaussian_inf_engine(bnet);

evidence = cell(1,N);
ym = [1 ; 1];
evidence{y} = Cy * ym;
[engine, ll] = enter_evidence(engine{1}, evidence);

marg = marginal_nodes(engine, x);
Exy = marg.mu;
Sxy = marg.Sigma;

% It results: 
% Exy = mx - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cy*ym - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cx*mx
% Sxy = Sx - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cx*Sx 

%%%%%%%%%%%%%%%%%%
% SECOND EXAMPLE %
%%%%%%%%%%%%%%%%%%
clear all
N = 2;
x = 1;  y = 2;
dag = zeros(N,N);
dag(1,2) = 1;

ns = [2 2];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% Cx x + Cy y = w
%
% x ~ N(mx, Sx)
% w ~ N(0 , R)

mx = [ 1 -1 ]';
Sx = [ 1 -1/2; -1/2 1];
R  = [ 1 0; 0 2];

Cx = [ 1 2  ; 4 5  ];
Cy = [ 1 0 1; 2 1 1];

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sx);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', [0 0]', 'cov', R, 'weights', -Cx);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
ym = [1 ; 1; -1];
evidence{y} = Cy * ym;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
Exy = marg.mu;
Sxy = marg.Sigma;

% It results: 
% mx - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cy*ym - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cx*mx
% Sx - Sx*Cx'*(R+Cx*Sx*Cx')^(-1)*Cx*Sx

%%%%%%%%%%%%%%%%%%
% THIRD  EXAMPLE %
%%%%%%%%%%%%%%%%%%

clear all

N = 2;
x = 1;  y = 2;
dag = zeros(N,N);
dag(1,2) = 1;

ns = [2 2];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% Cx x + Cy y = w
%
% x ~ N(mx, Sx)
% w ~ N(0 , R)

mx  = [ 0 0 ]';
R   = [ 2 1; 1 2];
Sx  = [ 100 0; 0 200];

Cx = [ 1 2  ; 4 5  ];
Cy = [ 0; 2];
ym  = 2;

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sx);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', [0 0]', 'cov', R, 'weights', Cx);

engine = gaussian_inf_engine(bnet);

evidence = cell(1,N);
evidence{y} = Cy * ym;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, x);
Exy = marg.mu;
Sxy = marg.Sigma;

% It results: 
Exy
(Sx^(-1)+Cx'*R^(-1)*Cx)^(-1)*Cx'*R^(-1)*Cy*ym
Sxy
(Sx^(-1)+Cx'*R^(-1)*Cx)^(-1)

%%%%%%%%%%%%%%%%%%%
% FOURTH  EXAMPLE %
%%%%%%%%%%%%%%%%%%%

clear all

N = 2;
x = 2;  y = 1;
dag = zeros(N,N);
dag(y,x) = 1;

ns = [1 2];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%  I x - Axy y = v
%
% x|y ~ N( Axy y , Sv)
%   y ~ N( my    , Sy)

mx  = [ 0 0 ]';
my  = [ 0 ]';
Sv  = [ 2 1; 1 2];
Sy  = eye(1);

Axy = [ 1 2 ]';

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sv, 'weights', Axy);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', my, 'cov', Sy);

engine = jtree_inf_engine(bnet);

evidence = cell(1,N);
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [y x]);
Exy = marg.mu;
Sxy = marg.Sigma;

% It results: 
% Sxy
% S = [Sy, Sy*Axy'; Axy*Sy, Sv + Axy*Sy*Axy'];
% S == Sxy

%%%%%%%%%%%%%%%%%%%
% FIFTH  EXAMPLE  %
%%%%%%%%%%%%%%%%%%%

clear all

N = 3;
x = 2;  y = 1;  z = 3;
dag = zeros(N,N);
dag(y,x) = 1;
dag(y,z) = 1;
dag(x,z) = 1;

ns = [2 2 2];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%  I x - Axy y         = v
%  I z - Bzx x - Bzy y = w
%
% x|y ~ N( Axy y , Sv)
%
% z|y ~ N( [Bzy Bzx] [y; x] , Sw)

mx  = [ 0 0 ]';
mz  = [ 0 0 ]';
my  = [ 0 0 ]';
Sv  = [ 2 1; 1 2];
Sw  = [ 100 0; 0 200];
Sy  = eye(2);

Axy = [ 1 2  ; 4 5  ];
Bzx = [ 0 2  ; 1 2  ];
Bzy = [ 0 2  ; 1 1  ];
zm  = [ 2 1 ]';

bnet.CPD{x} = gaussian_CPD(bnet, x, 'mean', mx, 'cov', Sv, 'weights', Axy);
bnet.CPD{z} = gaussian_CPD(bnet, z, 'mean', mz, 'cov', Sw, 'weights', [Bzy Bzx]);
bnet.CPD{y} = gaussian_CPD(bnet, y, 'mean', my, 'cov', Sy);

engine = jtree_inf_engine(bnet);

evidence     = cell(1,N);
evidence{z}  = zm;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [x y]);
Exy = marg.mu;
Sxy = marg.Sigma;

S = [Sy, Sy*Axy'; Axy*Sy, Sv + Axy*Sy*Axy'];
% It results: 
% Exy
% (Sx^(-1)+Cx'*R^(-1)*Cx)^(-1)*Cx'*R^(-1)*Cy*ym
% Sxy
% (Sx^(-1)+Cx'*R^(-1)*Cx)^(-1)

(S^(-1)+[Bzy Bzx]'*Sw^(-1)*[Bzy Bzx])^(-1)

%%%%%%%%%%%%%%%%%%%
% SIXTH  EXAMPLE  %
%%%%%%%%%%%%%%%%%%%

clear all

N = 3;
dx = 2;  dy = 1;  dz = 3;
dag = zeros(N,N);
dag(dy,dx) = 1;
dag(dy,dz) = 1;
dag(dx,dz) = 1;

ns = [2 2 2];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%  I dx + Dy dy         = v
%  I dz - Yx dx - Yy dy = w
%
%       dy ~ N(         my , Sy)
%
%    dx|dy ~ N( mx + Dy dy , Sv)
%
% dz|dx,dy ~ N( [Yy Yx] [dy; dx] , Sw)

mx  = [ 1 2 ]';
mz  = [ 3 4 ]';
my  = [ 5 6 ]';
Sv  = [ 2 1; 1 2];
Sw  = [ 100 0; 0 200];
Sy  = eye(2);

Dy = [ 1 2  ; 4 5  ];
Yx = [  0 2  ;  1 2  ];
Yy = [ -1 2  ;  1 1  ];
zm  = [ 2 1 ]';

bnet.CPD{dx} = gaussian_CPD(bnet, dx, 'mean', mx, 'cov', Sv, 'weights', -Dy);
bnet.CPD{dz} = gaussian_CPD(bnet, dz, 'mean', mz, 'cov', Sw, 'weights', [Yy Yx]);
bnet.CPD{dy} = gaussian_CPD(bnet, dy, 'mean', my, 'cov', Sy);

engine = jtree_inf_engine(bnet);

evidence     = cell(1,N);
evidence{dz}  = zm;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [dy dx]);
Exy_z = marg.mu;
Sxy_z = marg.Sigma;

% I dx + Dy dy = v
%
%    dy ~ N(         my , Sy)
% dx|dy ~ N( mx + Dy dy , Sv)
%
% Induces a distribution p(x,y)
%
% x,y ~ N(mxy, Sxy)
%
Sxy = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
mxy = -Sxy*[-Dy'*Sv^(-1)*mx - Sy^(-1)*my; -Sv^(-1)*mx];

% It results: 
Y = [Yy Yx];

Sxy_z
(Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)
Exy_z
mxy + (Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)*Y'*Sw^(-1)*((zm-mz)-Y*mxy) 



