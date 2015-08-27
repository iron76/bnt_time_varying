function CPD = update_ess_modified(CPD, CPDp, fmarginal, evidence, ns, cnodes, hidden_bitv)
% UPDATE_ESS Update the Expected Sufficient Statistics of a Gaussian node modified to support 
%            inertial paramters estimation 
% function CPD = update_ess(CPD, fmarginal, evidence, ns, cnodes, hidden_bitv)

%if nargin < 6
%  hidden_bitv = zeros(1, max(fmarginal.domain));
%  hidden_bitv(find(isempty(evidence)))=1;
%end

dom = fmarginal.domain;
self = dom(end);
ps = dom(1:end-1);
cps = myintersect(ps, cnodes);
dps = mysetdiff(ps, cps);

CPD.nsamples = CPD.nsamples + 1;            

% this is hardcoded for this node
ss = 6;
cpsz = 6;
dpsz = 1;

if (~any(any(CPD.WXsum)) && ~isempty(CPD.WXsum))
    CPD.WXsum  = zeros(ns(self),1);
end
if (~any(any(CPD.WXXsum)) && ~isempty(CPD.WXXsum))
    CPD.WXXsum = zeros(ns(self),ns(self));
end
if (~any(any(CPD.WXYsum)) && ~isempty(CPD.WXYsum))
    CPD.WXYsum = zeros(ns(self),ns(self));
end

if (~any(any(CPDp.WXsum)) && ~isempty(CPDp.WXsum))
    CPDp.WXsum  = zeros(ns(self),1);
end
if (~any(any(CPDp.WXXsum)) && ~isempty(CPDp.WXXsum))
    CPDp.WXXsum = zeros(ns(self),ns(self));
end
if (~any(any(CPDp.WXYsum)) && ~isempty(CPDp.WXYsum))
    CPDp.WXYsum = zeros(ns(self),ns(self));
end

% create data structure for inertial parameters ess
if ~(any(any(CPD.Asum)) && ~isempty(CPD.Asum))
    %disp('resetting CPD.Asum');
    CPD.Asum  = zeros(10,10);
end
if ~(any(any(CPD.Bsum)) && ~isempty(CPD.Bsum))
    %disp('resetting CPD.Bsum');
    CPD.Bsum = zeros(10,1);
end

if (~any(any(CPDp.Asum)) && ~isempty(CPDp.Asum))
    %disp('resetting CPDp.Asum');
    CPDp.Asum  = zeros(10,10);
end
if (~any(any(CPDp.Bsum)) && ~isempty(CPDp.Bsum))
    %disp('resetting CPDp.Bsum');
    CPDp.Bsum = zeros(10,1);
end

% Set weights and mean from node structure
% In this node we don't explicitly write the weights relating the 
% parent to the node and the mean of the gaussian CPD, but 
% we compute them from the given data (inertial parameters and twist)

% The weights with respect to the parent (spatial acceleration)
% are just the inertia matrix in link frame, that can be obtained from 
% inertial parameters
% B was CPD.weigths
% TODO \TODO check effect of weights  
B = inertiaMatrixFromInertialParams(CPD.inertial_params);

% currentMean was CPD.mean 
currentMean = crf(CPD.twist)*B*CPD.twist;

% Let X be the cts parent (if any), Y be the cts child (self).

if ~hidden_bitv(self) & ~any(hidden_bitv(cps)) & all(hidden_bitv(dps))
  % Speedup for the common case that all cts nodes are observed, all discrete nodes are hidden
  % Since X and Y are observed, SYY = 0, SXX = 0, SXY = 0
  % Since discrete parents are hidden, we do not need to add evidence to w.
  
  % we do not support this optimization for now in this type of node
  assert(false)
  
  w = fmarginal.T(:);
  CPD.Wsum = CPDp.Wsum + w;
  y = evidence{self} - currentMean;
  Cyy = y*y';
  if ~CPD.useC
     WY = repmat(w(:)',ss,1); % WY(y,i) = w(i)
     WYY = repmat(reshape(WY, [ss 1 dpsz]), [1 ss 1]); % WYY(y,y',i) = w(i)
     %CPD.WYsum = CPD.WYsum +  WY .* repmat(y(:), 1, dpsz);
     CPD.WYsum = CPDp.WYsum +  y(:) * w(:)';
     CPD.WYYsum = CPDp.WYYsum + WYY  .* repmat(reshape(Cyy, [ss ss 1]), [1 1 dpsz]);
  else
     W = w(:)';
     W2 = reshape(W, [1 1 dpsz]);
     CPD.WYsum = CPDp.WYsum +  rep_mult(W, y(:), size(CPD.WYsum)); 
     CPD.WYYsum = CPDp.WYYsum + rep_mult(W2, Cyy, size(CPD.WYYsum));
  end
  if cpsz > 0 % X exists

    x = cat(1, evidence{cps});
    z = B*x;
    Cxx = z*z';
    Cxy = z*y';
    WX = repmat(w(:)',ss,1); % WX(x,i) = w(i)
    WXX = repmat(reshape(WX, [ss 1 dpsz]), [1 ss 1]); % WXX(x,x',i) = w(i)
    WXY = repmat(reshape(WX, [ss 1 dpsz]), [1 ss 1]); % WXY(x,y,i) = w(i)
    if ~CPD.useC
      CPD.WXsum = CPDp.WXsum + WX .* repmat(z(:), 1, dpsz);
      CPD.WXXsum = CPDp.WXXsum + WXX .* repmat(reshape(Cxx, [ss ss 1]), [1 1 dpsz]);
      CPD.WXYsum = CPDp.WXYsum + WXY .* repmat(reshape(Cxy, [ss ss 1]), [1 1 dpsz]);
    else
      CPD.WXsum = CPDp.WXsum + rep_mult(W, z(:), size(CPD.WXsum));
      CPD.WXXsum = CPDp.WXXsum + rep_mult(W2, Cxx, size(CPD.WXXsum));
      CPD.WXYsum = CPDp.WXYsum + rep_mult(W2, Cxy, size(CPD.WXYsum));
    end
  end
  return;
end

% general (non-vectorized) case
fullm = add_evidence_to_gmarginal(fmarginal, evidence, ns, cnodes); % slow!

if dpsz == 1 % no discrete parents
  w = 1;
else
  w = fullm.T(:);
end

% ESS for covariance estimation
CPD.Wsum = CPDp.Wsum + w;
xi = 1:cpsz;
yi = (cpsz+1):(cpsz+ss);
for i=1:dpsz
  muY = fullm.mu(yi, i) - currentMean;
  SYY = fullm.Sigma(yi, yi, i);
  
  CPD.WYsum(:,i) = CPDp.WYsum(:,i) + w(i)*muY;
  CPD.WYYsum(:,:,i) = CPDp.WYYsum(:,:,i) + w(i)*(SYY + muY*muY'); % E[X Y] = Cov[X,Y] + E[X] E[Y]
  if cpsz > 0
    muZ = B*fullm.mu(xi, i);
    SZZ = B*fullm.Sigma(xi, xi, i)*B';
    SZY = B*fullm.Sigma(xi, yi, i);
    CPD.WXsum(:,i) = CPDp.WXsum(:,i) + w(i)*muZ;
    CPD.WXXsum(:,:,i) = CPDp.WXXsum(:,:,i) + w(i)*(SZZ + muZ*muZ');
    CPD.WXYsum(:,:,i) = CPDp.WXYsum(:,:,i) + w(i)*(SZY + muZ*muY');
  end
end  

% ESS for inertial parameters estimation
SAA = fullm.Sigma(xi, xi, 1);
muA = fullm.mu(xi, 1);
muAnotWeighted = diag(1./CPD.acceleration_weights)*muA;
EAA = SAA + muA*muA';
wrench_weights_diag = diag(CPD.wrench_weights);

invSigma_fiB = inv(CPD.cov);
% we can move this outside the samples loop to improve perforamnce
AquadForm = getQuadraticFormsForInertialParamsAESS(invSigma_fiB,wrench_weights_diag);
AquadFormSum = zeros(10,10);
for row =1:10
    for col = 1:10
        AquadFormSum(row,col) = trace(AquadForm{row,col}*EAA);
    end
end

SFA = fullm.Sigma(yi, xi, 1);
muF = fullm.mu(yi, 1);
EFA = SFA + muF*muA';

% we can move this outside the samples loop to improve perforamnce
BquadForm = getQuadraticFormsForInertialParamsBESS(invSigma_fiB,wrench_weights_diag);
BquadFormSum = zeros(10,1);
for i = 1:10
    BquadFormSum(i) = trace(BquadForm{i}*EFA);
end


coriolisRegressor = crf(CPD.twist)*inertiaRegressor(CPD.twist);
coriolisRegressorWeighted = wrench_weights_diag*coriolisRegressor;
inertiaRegressorAccWeighted =  wrench_weights_diag*inertiaRegressor(muA);

CPD.Asum = CPDp.Asum + AquadFormSum  + (coriolisRegressorWeighted)'*invSigma_fiB*inertiaRegressorAccWeighted + ...
          inertiaRegressorAccWeighted'*invSigma_fiB*coriolisRegressorWeighted + coriolisRegressorWeighted'*invSigma_fiB*coriolisRegressorWeighted;
%fprintf('CPDp.Bsum\n')
%display(CPDp.Bsum)
CPD.Bsum = CPDp.Bsum + BquadFormSum  + (coriolisRegressorWeighted)'*invSigma_fiB*muF;
%fprintf('BquadFormSum\n')
%display(BquadFormSum)
%fprintf('Other parts\n')
%display((coriolisRegressorWeighted)'*invSigma_fiB*muF)

