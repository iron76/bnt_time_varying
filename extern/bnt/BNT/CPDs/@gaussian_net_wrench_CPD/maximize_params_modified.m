function CPD = maximize_params_modified(CPD, temp)
% MAXIMIZE_PARAMS Set the params of a CPD to their ML values (Gaussian)
% CPD = maximize_params(CPD, temperature)
%
% Temperature is currently ignored.

if ~adjustable_CPD(CPD), return; end

cl_mean = zeros(6,1);
cl_weights = zeros(6,6);

if CPD.clamped_cov
  cl_cov = CPD.cov;
else
  cl_cov = [];
end


% [ssz psz Q] = size(CPD.weights);

% [ss cpsz dpsz] = size(CPD.weights); % ss = self size = ssz

ss = 6;
ssz = 6;
cpsz = 6;
psz = 6;
Q = 1;
dpsz = 1;

prior =  repmat(CPD.cov_prior_weight*eye(ssz,ssz), [1 1 Q]);

% Estimate covariance
[mean, CPD.cov, weights] = ...
    clg_Mstep_modified(CPD.Wsum, CPD.WYsum, CPD.WYYsum, [], CPD.WXsum, CPD.WXXsum, CPD.WXYsum, ...
	      'cov_type', CPD.cov_type, 'clamped_mean', cl_mean, ...
	      'clamped_cov', cl_cov, 'clamped_weights', cl_weights, ...
	      'tied_cov', CPD.tied_cov, ...
	      'cov_prior', prior);
      
% Estimate inertial parameters
CPD.inertial_params = CPD.Asum\CPD.Bsum;

% Bug fix 11 May 2003 KPM
% clg_Mstep collapses all discrete parents into one mega-node
% but convert_to_CPT needs access to each parent separately
sz = CPD.sizes;
ss = sz(end);

% Bug fix KPM 20 May 2003: 
cpsz = sum(sz(CPD.cps));
%if isempty(CPD.cps)
%  cpsz = 0;
%else
%  cpsz = sz(CPD.cps);
%end
dpsz = sz(CPD.dps);
CPD.cov = myreshape(CPD.cov, [ss ss dpsz]);
