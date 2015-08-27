function CPD = reset_ess(CPD)
% RESET_ESS Reset the Expected Sufficient Statistics for a Gaussian CPD.
% CPD = reset_ess(CPD)

CPD.nsamples = 0;    
CPD.Wsum = zeros(size(CPD.Wsum));
CPD.WYsum = zeros(size(CPD.WYsum));
CPD.WYYsum = zeros(size(CPD.WYYsum));
CPD.WXsum = zeros(size(CPD.WXsum));
CPD.WXXsum = zeros(size(CPD.WXXsum));
CPD.WXYsum = zeros(size(CPD.WXYsum));

CPD.Asum  = zeros(10,10);
CPD.Bsum  = zeros(10,1);
