function [m, C, W] = gaussian_net_wrench_get_canonical_parameters(CPD)
%GAUSSIANNETWRENCHGETCANONICALPARAMETERS Get the mean, weight and
%covariance of this CPD, encapsulating its internal representation
wrench_weights_diag = diag(CPD.wrench_weights);
inv_acc_weights_diag = diag(1./CPD.acceleration_weights);

B = inertiaMatrixFromInertialParams(CPD.inertial_params);
W = wrench_weights_diag*B*inv_acc_weights_diag;
C = CPD.cov(:,:);
m = wrench_weights_diag*crf(CPD.twist)*B*CPD.twist;

end

