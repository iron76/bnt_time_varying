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

%%%%%%%%%%%%%%%%%%%%%
% SEVENTS  EXAMPLE  %
%%%%%%%%%%%%%%%%%%%%%

clear all

N = 3;
dx1 = 2;  dx2 = 1;  dy = 3;
dag = zeros(N,N);
dag(dx2,dx1) = 1;
dag(dx2,dy)  = 1;
dag(dx1,dy)  = 1;

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

mx1  = [ 1 2 ]';
my   = [ 3 4 ]';
mx2  = [ 5 6 ]';
Sx1  = [ 2 1; 1 2];
Sy   = [ 100 0; 0 200];
Sx2  = eye(2);

D2 = [ 1 2  ; 4 5  ];
Y1 = [  0 2  ;  1 2  ];
Y2 = [ -1 2  ;  1 1  ];
ym  = [ 2 1 ]';

bnet.CPD{dx1} = gaussian_CPD(bnet, dx1, 'mean', mx1, 'cov', Sx1, 'weights', -D2);
bnet.CPD{dy } = gaussian_CPD(bnet, dy,  'mean', my,  'cov', Sy,  'weights', [Y2 Y1]);
bnet.CPD{dx2} = gaussian_CPD(bnet, dx2, 'mean', mx2, 'cov', Sx2);

engine = jtree_inf_engine(bnet);

evidence     = cell(1,N);
evidence{dy} = ym;
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [dx2 dx1]);
Ex_y = marg.mu;
Sx_y = marg.Sigma;

% I dx + Dy dy = v
%
%    dy ~ N(         my , Sy)
% dx|dy ~ N( mx + Dy dy , Sv)
%
% Induces a distribution p(x,y)
%
% x,y ~ N(mxy, Sxy)
%
% Sx  = [Sx2, -Sx2*D2'; -D2*Sx2, Sx1 + D2*Sx2*D2'];
Sx  = [Sx1 + D2*Sx2*D2', -D2*Sx2; -Sx2*D2', Sx2];
% mx  = -Sx*[-D2'*Sx1^(-1)*mx1 - Sx2^(-1)*mx2; -Sx1^(-1)*mx1];
mx  = -Sx*[-Sx1^(-1)*mx1; -D2'*Sx1^(-1)*mx1 - Sx2^(-1)*mx2];

% It results: 
Y = [Y1 Y2];

mx
[mx1-D2*mx2; mx2]
[Sx_y(3:4, 3:4) Sx_y(3:4, 1:2); Sx_y(1:2, 3:4) Sx_y(1:2, 1:2)]
(Sx^(-1)+Y'*Sy^(-1)*Y)^(-1)
[Ex_y(3:4,1); Ex_y(1:2,1)]
mx + (Sx^(-1)+Y'*Sy^(-1)*Y)^(-1)*Y'*Sy^(-1)*((ym-my)-Y*mx) 

%%%%%%%%%%%%%%%%%%%
% EIGTH  EXAMPLE  %
%%%%%%%%%%%%%%%%%%%

clear all

N = 3;
 d = 1;  e = 2;   y = 3;
sd = 4; se = 3;  sy = 2;
dag = zeros(N,N);
dag(d,e) = 1;
dag(d,y) = 1;

ns = [sd se sy];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%  We consider the case:
%                  d ~ N(md, Sd)
%  D d + bD      = v ~ N(me, Se)
%  Y d + bY - y  = w ~ N(my, Sy)
%
%  The analytic solution is:
%  d|y ~ (mu, S)
%  S  = (Sd_inv + D' SD_inv D + Y' SY_inv Y)^-1
%  mu = S(Y' SY_inv (y-bY) + Sd_inv md + D' SD_inv (mD - bD))

md  = [ 1 2 3 4]';
me  = [ 1 2 3 ]';
my  = [ 5 6 ]';
Sd  = blkdiag([ 2 1; 1 2], [ 1 0; 1 2]);
Se  = [ 100 10 10 ; 10 200 0; 10 0 20];
Sy  = eye(2);

D   = [ 1 2  0 1 ; 4 5 1 0; -1 0 0 1];
Y   = [ 0 2  0 1 ; 1 2 1 0];

bD  = [ 1; -1; 0];
bY  = [-1; -2];

ym  = [1; 1];

bnet.CPD{d}  = gaussian_CPD(bnet, d, 'mean',    md, 'cov', Sd);
bnet.CPD{e}  = gaussian_CPD(bnet, e, 'mean', me+bD, 'cov', Se, 'weights', D);
bnet.CPD{y}  = gaussian_CPD(bnet, y, 'mean', my+bY, 'cov', Sy, 'weights', Y);


engine = jtree_inf_engine(bnet);

evidence     = cell(1,N);
evidence{y}  = ym;
evidence{e}  = zeros(se, 1);


[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, d);
Ed_y = marg.mu;
Sd_y = marg.Sigma;


Sd_y
(Sd^(-1) + D' * Se^(-1) * D + Y' * Sy^(-1) * Y)^(-1)
Ed_y
Sd_y * (Y' * Sy^(-1) * (ym-my-bY) + Sd^(-1) * md + D'* Se^(-1) * (- me - bD))

%%%%%%%%%%%%%%%%%%%
% NINETH  EXAMPLE %
%%%%%%%%%%%%%%%%%%%

clear all

N = 3;
 dy = 1;  dx = 2;   y = 3;
sdy = 1; sdx = 3;  sy = 2;
dag = zeros(N,N);
dag(dy,dx) = 1;
dag(dy,y)  = 1;
dag(dx,y)  = 1;

ns = [sdy sdx sy];   % vector-valued

dnodes = [];  % no discrete nodes
bnet = mk_bnet(dag, ns, 'discrete', dnodes);

%  We consider the case:
%
%  d = [dx dy]              dy ~ N(mdy, Sdy)
%  Dx dx + Dy dy + bD      = v ~ N(me, Se)
%  Yx dx + Yy dy + bY - y  = w ~ N(my, Sy)
%
%  The analytic solution is:
%  d|y ~ (mu, S)
%  S  = (Sd_inv + D' SD_inv D + Y' SY_inv Y)^-1
%  mu = S(Y' SY_inv (y-bY) + Sd_inv md + D' SD_inv (mD - bD))
%
%  In a Bayesian network, it can be represented as follows:
%
%  d = [dx dy]                          dy ~ N(mdy, Sdy)
%     dx + Dx^-1 Dy dy + Dx^-1 bD      = v ~ N(Dx^-1 me, Dx^-1 Se Dx^-T)
%  Yx dx +       Yy dy +       bY - y  = w ~ N(my, Sy)
%

mdy = 4;
me  = [ 1 2 3 ]';
my  = [ 5 6 ]';
Sdy = 2;
Se  = [ 100 10 10 ; 10 200 0; 10 0 20];
Sy  = eye(2);

Sd_inv  = blkdiag(zeros(sdx,sdx), inv(Sdy));
md      = [zeros(sdx,1); mdy];

D   = [ 1 2  0 1 ; 4 5 1 0; -1 0 0 1];
Y   = [ 0 2  0 1 ; 1 2 1 0];

Dx  = D(:,     1:sdx);
Dy  = D(:, sdx+1:end);
Yx  = Y(:,     1:sdx);
Yy  = Y(:, sdx+1:end);

bD  = [ 1; -1; 0];
bY  = [-1; -2];

ym  = [1; 1];

bnet.CPD{dy}  = gaussian_CPD(bnet, dy, 'mean',   mdy, 'cov', Sdy);
bnet.CPD{dx}  = gaussian_CPD(bnet, dx, 'mean', inv(Dx) * (-me-bD), 'cov', inv(Dx) * Se * inv(Dx'), 'weights', - inv(Dx) * Dy);
bnet.CPD{y}   = gaussian_CPD(bnet,  y, 'mean', my+bY, 'cov', Sy, 'weights', [Yy Yx]);


engine = jtree_inf_engine(bnet);

evidence     = cell(1,N);
evidence{y}  = ym;

[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, [dy dx]);
Ed_y = marg.mu;
Sd_y = marg.Sigma;


Sd_y([sdy+(1:sdx) 1:sdy], [sdy+(1:sdx) 1:sdy])
(Sd_inv + D' * Se^(-1) * D + Y' * Sy^(-1) * Y)^(-1)
Ed_y([sdy+(1:sdx) 1:sdy])
(Sd_inv + D' * Se^(-1) * D + Y' * Sy^(-1) * Y)^(-1) * (Y' * Sy^(-1) * (ym-my-bY) + Sd_inv * md + D'* Se^(-1) * (- me - bD))

