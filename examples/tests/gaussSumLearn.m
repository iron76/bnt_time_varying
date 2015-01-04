clc
close
clear

N = 3;
dx = 2;  dy = 1;  dz = 3;
dag = zeros(N,N);
dag(dy,dx) = 1;
dag(dy,dz) = 1;
dag(dx,dz) = 1;

ns = [2 5 3];   % vector-valued

dnodes    = [];         % no discrete nodes
obs_nodes = [dz];    % measured nodes
bnet      = mk_bnet(dag, ns, 'discrete', dnodes, 'observed', obs_nodes);

%  I dx + Dy dy         = v
%  I dz - Yx dx - Yy dy = w
%
%    dy    ~ N(         my , Sy)
%
% dx|dy    ~ N( mx + Dy dy , Sv)
%
% dz|dx,dy ~ N( [Yy Yx] [dy; dx] , Sw)

mx  = [ 1 2 3 4 5 ]';
mz  = [ 3 4 1 ]';
my  = [ 5 6 ]';
Sv  = eye(5)*1e-5;
Sw  = eye(3);
Sy  = eye(2);

Dy = [ 1 2  ; 4 5 ; 1 1; 0 1; 1 0 ];
Yx = [  0 2 1 -1 0 ;  1 2 1 -1 0; 1 2 1 4 1];
Yy = [ -1 2  ;  1 1 ; 0 1 ];
zm  = [ 2 1 0 ]';

bnet.CPD{dx} = gaussian_CPD(bnet, dx, 'mean', mx, 'cov', Sv, 'weights', -Dy, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
bnet.CPD{dz} = gaussian_CPD(bnet, dz, 'mean', mz, 'cov', Sw, 'weights', [Yy Yx], 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 1e-12);
bnet.CPD{dy} = gaussian_CPD(bnet, dy, 'mean', my, 'cov', Sy, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);

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
% with:
%
% Sxy = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
% mxy = -Sxy*[-Dy'*Sv^(-1)*mx - Sy^(-1)*my; -Sv^(-1)*mx];
% 
% Given the following:
%
%  I dz - Yx dx - Yy dy = w
%
% dz|dx,dy ~ N( [Yy Yx] [dy; dx] , Sw)
%
% gives the distribution p(x,y|z)
%
% x,y|z ~ N(mxy_z, Sxy_z)
%
% with:
%
% Y     = [Yy Yx];
% Sxy_z = (Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)
% mxy_z = mxy + (Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)*Y'*Sw^(-1)*((zm-mz)-Y*mxy)

% Numerical verification:
%
% Y = [Yy Yx];
% Sxy_z
% Sxy = [Sy, -Sy*Dy'; -Dy*Sy, Sv + Dy*Sy*Dy'];
% (Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)
% Exy_z
% mxy = -Sxy*[-Dy'*Sv^(-1)*mx - Sy^(-1)*my; -Sv^(-1)*mx];
% mxy + (Sxy^(-1)+Y'*Sw^(-1)*Y)^(-1)*Y'*Sw^(-1)*((zm-mz)-Y*mxy) 

% Learning stage
n        = 1000;
samples  = cell(N,n);
evidence = cell(1,N);

% Perturb the initial estimation
bnet0          = bnet;
bnet0.CPD{dx}  = set_fields(bnet.CPD{dx}, 'cov',     get_field(bnet.CPD{dx}, 'cov'));
bnet0.CPD{dz}  = set_fields(bnet.CPD{dz}, 'cov',     get_field(bnet.CPD{dz}, 'cov'));
bnet0.CPD{dy}  = set_fields(bnet.CPD{dy}, 'cov',     get_field(bnet.CPD{dy}, 'cov'));
engine         = jtree_inf_engine(bnet0);
engine         = enter_evidence(engine, evidence); 

i_hidden       = setdiff(1:N, bnet.observed);
for i = 1 : n
  samples(:,i) = sample_bnet(bnet, 'evidence', evidence);
  for j = 1 : length(i_hidden)
     samples{i_hidden(j),i} = [];
  end
end

bnet_hat = learn_params_em(engine, samples, 90, 1e-5);

% Visualize results 
err0 = norm(get_field(bnet.CPD{dz}, 'cov') - get_field(bnet0.CPD{dz}, 'cov'));
errH = norm(get_field(bnet.CPD{dz}, 'cov') - get_field(bnet_hat.CPD{dz}, 'cov'));
if (err0 > errH)
   fprintf('[INF] The variance estimation improved from %f to %f. \n', err0, errH)
else
   fprintf('[ERR] The variance estimation did not improved! \n')
end
