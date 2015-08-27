function res = testRegressorsExpectations

tol = 1e-5;

% Test for checking that the formula used to compute the expectation 
% of the quantities involving regressors in the M step of the modified 
% EM for inertial parameters estimation
res = 0;

% to check the computations, we will compare the expectation formulas with
% the case where a (the parent distribution variable) is perfectly known,
% and then E(a*a^\top) is simply given as a*a^\top
% create random spatial acceleration and velocity
a = rand(6,1);
v = rand(6,1);
fiB = rand(6,1);
wrench_weights_diag = diag(rand(6,1));
acc_weights_diag    = diag(rand(6,1));
a_weighted = acc_weights_diag*a;
fiB_weighted = wrench_weights_diag*fiB;
invSigma_fiB = diag(rand(6,1));
EAA= a_weighted*a_weighted';

coriolisRegressor = crf(v)*inertiaRegressor(v);
coriolisRegressorWeighted = wrench_weights_diag*coriolisRegressor;
inertiaRegressorAccWeighted =  wrench_weights_diag*inertiaRegressor(a_weighted);

% we first compare the Y^\top Sigma^-1 Y term
completeRegressorWeighted = inertiaRegressorAccWeighted + coriolisRegressorWeighted;
YTinvSigmaYbaseLine = completeRegressorWeighted'*invSigma_fiB*completeRegressorWeighted;

AquadForm = getQuadraticFormsForInertialParamsAESS(invSigma_fiB,wrench_weights_diag);
AquadFormSum = zeros(10,10);
for row =1:10
    for col = 1:10
        % Check https://en.wikipedia.org/wiki/Quadratic_form_%28statistics%29
        AquadFormSum(row,col) = trace(AquadForm{row,col}*EAA);
    end
end

if( not(all(isalmost(AquadFormSum,inertiaRegressorAccWeighted'*invSigma_fiB*inertiaRegressorAccWeighted,tol))) ) 
    disp('Something wrong in getQuadraticFormsForInertialParamsAESS');
    res = 1;
    return
end


YTinvSigmaY = AquadFormSum  + (coriolisRegressorWeighted)'*invSigma_fiB*inertiaRegressorAccWeighted + ...
          inertiaRegressorAccWeighted'*invSigma_fiB*coriolisRegressorWeighted + coriolisRegressorWeighted'*invSigma_fiB*coriolisRegressorWeighted;

if( not(all(isalmost(YTinvSigmaY,YTinvSigmaYbaseLine,tol))) ) 
    disp('Something wrong in getQuadraticFormsForInertialParamsAESS');
    res = 1;
end    

% we then compare the  Y^\top Sigma^-1 fBi term
YTinvSigmafiBbaseLine = completeRegressorWeighted'*invSigma_fiB*fiB_weighted;

EFA = fiB_weighted*a_weighted';

BquadForm = getQuadraticFormsForInertialParamsBESS(invSigma_fiB,wrench_weights_diag);
BquadFormSum = zeros(10,1);
for i = 1:10
    % Check https://en.wikipedia.org/wiki/Quadratic_form_%28statistics%29
    BquadFormSum(i) = trace(BquadForm{i}*EFA);
end

YTinvSigmafiB = BquadFormSum  + (coriolisRegressorWeighted)'*invSigma_fiB*fiB_weighted;

if( not(all(isalmost(BquadFormSum,inertiaRegressorAccWeighted'*invSigma_fiB*fiB_weighted,tol))) ) 
    disp('Something wrong in getQuadraticFormsForInertialParamsBESS');
    res = 1;
    return
end

if( not(all(isalmost(YTinvSigmafiB,YTinvSigmafiBbaseLine,tol))) ) 
    disp('Something wrong in getQuadraticFormsForInertialParamsBESS');
    res = 1;
end  

end