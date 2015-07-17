function CPD = gaussian_net_wrench_CPD(bnet, self, varargin)
% GAUSSIAN_NET_WRENCH_CPD Make a conditional linear Gaussian distrib,
%                         following the bayesian dependency between the spatial 
%                         acceleration of a link and its net wrench (i.e. the derivative 
%                         of the spatial momentum of the link)
%
% CPD = gaussian_CPD(bnet, node, ...) will create a CPD with random parameters,
% where node is the number of a node in this equivalence class.

% This CPD can only have a single continuous parent of 6 elements, 
% the spatial acceleration of the link expressed in link frame. 
% Given the spatial acceleration a, the gaussian distribution for this 
% CPD is given by: 
%  Y|X=a ~ N( I*a + v*I*v , Sigma)
% Where v is the spatial twist of the link and I is the spatial inertia. 
% We define this special CPD because during learning we want to learn the 
% inertial parameters in I exploiting its structure, i.e. given v we want to 
% learn both Sigma and the phi such that I*a+v*I*v = Regr(a,0)\phi + Regr(0,v)\phi 
%
% The list below gives optional arguments [default value in brackets].
% (Let ns(i) be the size of node i, ns(X) and ns(Y) are bounded to be 6 )
%
% twist      - twist(:) is the twist of the rigid body [ zeros(6,1) ]
% cov        - Sigma(:,:) is the covariance [ 100*eye(6,6) ]
% cov_type   - if 'diag', Sigma(:,:) is diagonal [ 'full' ]
% clamp_cov  - if 1, we do not adjust Sigma(:,:) during learning [0]
% cov_prior_weight - weight given to I prior for estimating Sigma [0.01]
% cov_prior_entropic - if 1, we also use an entropic prior for Sigma [0]
%
% e.g., CPD = gaussian_net_wrench_CPD(bnet, i, 'twist', [0.4,0.5,0.5,1.3,5.6,5.6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% STILL NEED TO MODIFY EVERYTHING AFTER THIS LINE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
  % This occurs if we are trying to load an object from a file.
  CPD = init_fields;
  clamp = 0;
  CPD = class(CPD, 'gaussian_net_wrench_CPD', generic_CPD(clamp));
  return;
elseif isa(bnet, 'gaussian_net_wrench_CPD')
  % This might occur if we are copying an object.
  CPD = bnet;
  return;
end
CPD = init_fields;
 
CPD = class(CPD, 'gaussian_net_wrench_CPD', generic_CPD(0));

args = varargin;
ns = bnet.node_sizes;
ps = parents(bnet.dag, self);
dps = myintersect(ps, bnet.dnodes);
cps = myintersect(ps, bnet.cnodes);
fam_sz = ns([ps self]);

CPD.self = self;
CPD.sizes = fam_sz;

% Figure out which (if any) of the parents are discrete, and which cts, and how big they are
% dps = discrete parents, cps = cts parents
CPD.cps = find_equiv_posns(cps, ps); % cts parent index
CPD.dps = find_equiv_posns(dps, ps);
ss = fam_sz(end);
psz = fam_sz(1:end-1);
dpsz = prod(psz(CPD.dps));
cpsz = sum(psz(CPD.cps));

% set default params
CPD.mean = randn(ss, dpsz);
CPD.cov = 100*repmat(eye(ss), [1 1 dpsz]);    
CPD.weights = randn(ss, cpsz, dpsz);
CPD.cov_type = 'full';
CPD.tied_cov = 0;
CPD.clamped_mean = 0;
CPD.clamped_cov = 0;
CPD.clamped_weights = 0;
CPD.cov_prior_weight = 0.01;
CPD.cov_prior_entropic = 0;
nargs = length(args);
if nargs > 0
  CPD = set_fields(CPD, args{:});
end

% Make sure the matrices have 1 dimension per discrete parent.
% Bug fix due to Xuejing Sun 3/6/01
CPD.mean = myreshape(CPD.mean, [ss ns(dps)]);
CPD.cov = myreshape(CPD.cov, [ss ss ns(dps)]);
CPD.weights = myreshape(CPD.weights, [ss cpsz ns(dps)]);

% Precompute indices into block structured  matrices
% to speed up CPD_to_lambda_msg and CPD_to_pi
cpsizes = CPD.sizes(CPD.cps);
CPD.cps_block_ndx = cell(1, length(cps));
for i=1:length(cps)
  CPD.cps_block_ndx{i} = block(i, cpsizes);
end

%%%%%%%%%%% 
% Learning stuff

% expected sufficient statistics 
CPD.Wsum = zeros(dpsz,1);
CPD.WYsum = zeros(ss, dpsz);
CPD.WXsum = zeros(cpsz, dpsz);
CPD.WYYsum = zeros(ss, ss, dpsz);
CPD.WXXsum = zeros(cpsz, cpsz, dpsz);
CPD.WXYsum = zeros(cpsz, ss, dpsz);

% For BIC
CPD.nsamples = 0;
switch CPD.cov_type
 case 'full',
  % since symmetric 
    %ncov_params = ss*(ss-1)/2; 
    ncov_params = ss*(ss+1)/2; 
  case 'diag',
    ncov_params = ss;
  otherwise
    error(['unrecognized cov_type ' cov_type]);
end
% params = weights + mean + cov
if CPD.tied_cov
  CPD.nparams = ss*cpsz*dpsz + ss*dpsz + ncov_params;
else
  CPD.nparams = ss*cpsz*dpsz + ss*dpsz + dpsz*ncov_params;
end

% for speeding up maximize_params
CPD.useC = exist('rep_mult');

clamped = CPD.clamped_mean & CPD.clamped_cov & CPD.clamped_weights;
CPD = set_clamped(CPD, clamped);

%%%%%%%%%%%

function CPD = init_fields()
% This ensures we define the fields in the same order 
% no matter whether we load an object from a file,
% or create it from scratch. (Matlab requires this.)

CPD.self = [];
CPD.sizes = [];
CPD.twist = [];
CPD.cov = [];
CPD.clamped_cov = [];
CPD.cov_type = [];
CPD.tied_cov = [];
CPD.Wsum = [];
CPD.WYsum = [];
CPD.WXsum = [];
CPD.WYYsum = [];
CPD.WXXsum = [];
CPD.WXYsum = [];
CPD.nsamples = [];
CPD.nparams = [];            
CPD.cov_prior_weight = [];
CPD.cov_prior_entropic = [];
CPD.useC = [];
CPD.cps_block_ndx = [];
