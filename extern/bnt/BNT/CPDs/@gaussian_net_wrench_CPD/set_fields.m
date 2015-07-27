function CPD = set_fields(CPD, varargin)
% SET_PARAMS Set the parameters (fields) for a gaussian_CPD_net_wrench object
% CPD = set_params(CPD, name/value pairs)
%
% The following optional arguments can be specified in the form of name/value pairs:
%
% cov        - Sigma(:,:) is the covariance
% cov_type   - if 'diag', Sigma(:,:,i) is diagonal 
% tied_cov   - if 1, we constrain Sigma(:,:,i) to be the same for all i
% clamp_cov  - if 1, we do not adjust Sigma(:,:,i) during learning 
% clamp      - if 1, we do not adjust any params
% cov_prior_weight - weight given to I prior for estimating Sigma
% cov_prior_entropic - if 1, we also use an entropic prior for Sigma [0]
% inertial_params = inertial_params(:) is the vector of inertial parameters
% clamp_inertial_params = if 1, we do not adjust inertial_params during learning 
% force_weights =  the vector of weights used for the wrench
% acceleration_weights = the vector of weights used for accelerations
%
% e.g., CPD = set_params(CPD, 'mean', [0;0])

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'twist',        CPD.twist = args{i+1}; 
   case 'cov',         CPD.cov = args{i+1}; 
   case 'cov_type',    CPD.cov_type = args{i+1}; 
   case 'tied_cov',    CPD.tied_cov = args{i+1};
   case 'clamp_cov',   CPD.clamped_cov = args{i+1};
   case 'clamp',  clamp = args{i+1};
    CPD.clamped_mean = clamp;
    CPD.clamped_cov = clamp;
    CPD.clamped_weights = clamp;
   case 'cov_prior_weight',  CPD.cov_prior_weight = args{i+1};
   case 'cov_prior_entropic',  CPD.cov_prior_entropic = args{i+1};
   case 'inertial_params', CPD.inertial_params = args{i+1};
   case 'clamp_inertial_params', CPD.clamped_inertial_params = args{i+1};
   case 'wrench_weights', CPD.wrench_weights = args{i+1};
   case 'force_weights', CPD.wrench_weights = args{i+1};
   case 'acceleration_weights', CPD.acceleration_weights = args{i+1}; 
   otherwise,  
    error(['invalid argument name ' args{i}]);
  end
end
