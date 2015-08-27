function y = sample_node(CPD, pev)
% SAMPLE_NODE Draw a random sample from P( f | a )  (gaussian)
% y = sample_node(CPD, parent_evidence)
%
% pev{i} is the value of the i'th parent (if there are any parents)
% y is the sampled value (a scalar or vector)

wrench_weights_diag = diag(CPD.wrench_weights);
inv_acc_weights_diag = diag(1./CPD.acceleration_weights);

isempty(CPD.dps);
assert( ~(isempty(CPD.cps)) );
  pev = pev(:);
  x = pev{1};
  B = inertiaMatrixFromInertialParams(CPD.inertial_params);
  % currentMean was CPD.mean 
  currentMean = wrench_weights_diag*(B*inv_acc_weights_diag*x+crf(CPD.twist)*B*CPD.twist);
  y = gsamp(currentMean, CPD.cov(:,:), 1);
  y = y(:);
end
