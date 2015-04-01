function res = testCalibration(dmodel_DNEA, ymodel_DNEA, dmodel, ymodel, S_dmodel)

% This function checks the calibration procedure. Given x, the calibration
% consists in estimating [d; dx] which satisy the following:
%
% D(x+dx) d + b(x+dx) = 0
% Y(x+dx) d - y    = 0
%
% The estimation is similar to a Newton-like optimization step, starting
% from an initial estimation [d_bar; dx=0]. The first order
% approximation of the above non-linear equations around [d_bar; dx=0]
% is given by:
%
% D(x) d + b(x) + dDb(d_bar, x) dx = 0
% Y(x) d - y    + dY(x)         dx = 0
%
% where dDb(d_bar, x) is the derivative of D(x+dx)d+b with respect to dx
% evaluated at [d_bar; dx=0] and dY(x_bar) is the derivative of Y(x+dx) with
% respect to x evaluated at dx=0. In order to check the procedure, we
% start by assuming [d_bar; dx=0] such that it satisfies the non-linear
% equation system. The first set of equations is guranteed by choosing:
%
% d_bar = myRNEA.d;
%
% The second set of equations is guaranteed by choosing:
%
% y = mySNEA.simY(d_bar);
%
% In the specific implmentation Y is not a function of x and therefore the
% linerarization simplifies as follows:
%
% D(x) d + b(x) + dDb(d_bar, x) dx = 0
% Y(x) d - y                       = 0
%
% The estimation for [d; dx] is obtained by solving the above equations in
% [d;dx], thus it consits of solving a linear system. If the matrix:
%
% | D(x)    dDb(d_bar, x) |
% | Y(x)    0             |
%
% is full-rank, a unique solution exists and should coincide with
% [d; dx] = [d_bar; 0].


ymodel_RNEA = autoSensRNEA(dmodel);
ymodel_RNEA = autoSensStochastic(ymodel_RNEA, 1e-5);

res = 0;
NB  = dmodel.NB;
T   = 20;

Q         = .1*eye(2*NB);
Sx        = diag(1e4*S_dmodel*rand(NB*2, 1));
eq        = 20*randn(NB,1)*pi/180;
edq       = 20*randn(NB,1)*pi/180;
dx        = zeros(2*NB,1);
q         = randn(NB,1)*0.5;
dq        = randn(NB,1)*0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_RNEA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_DNEA);
myDNEA  = DNEA(myModel, mySens);
% myModel = myModel.updateVarianceDNEA(0, myDNEA.dDb_s);
% mySens  = mySens.updateVarianceDNEA(0, myDNEA.dby_s);
% myDNEA  = DNEA(myModel, mySens);

for i = 1 : T
   myRNEA    = myRNEA.setState(q,dq);
   myRNEA    = myRNEA.solveID();
   d         = myRNEA.d;
   
   d_pri  = d + randn(size(d)).*abs(d).*0.1;
   d_bar  = d_pri;
   
   x      = [q; dq];
   x_bar  = [q; dq] + [eq; edq];
   x_pri  = x_bar + dx;
   
   disp(['Distance forom real x is: ' num2str(norm(x_pri - x))]);
      
   myDNEA = myDNEA.setState(x(1:NB ,1),x(NB+1:end,1));
   y      = myDNEA.simY(d,  x(1:NB ,1),x(NB+1:end,1));
   
   myDNEA = myDNEA.setState(x_pri(1:NB ,1),x_pri(NB+1:end,1));
   myDNEA = myDNEA.setY(y);
   
   myDNEA = myDNEA.setD(d_bar);
   myDNEA = myDNEA.setDprior(d_pri);
   myDNEA = myDNEA.setXprior(x_pri);
   myDNEA = myDNEA.setXvariance(Sx);
   myDNEA = myDNEA.solveID();
   
   errDNEAd(i)  = norm(myDNEA.d - d);
   errDNEAq(i)  = norm(myDNEA.x(   1:NB ,1)  -  q);
   errDNEAdq(i) = norm(myDNEA.x(NB+1:end,1)  - dq);
   
   errd(i)  = norm(d-d_bar);
   errq(i)  = norm(eq);
   errdq(i) = norm(edq);
   errdx(i) = norm(dx+[eq; edq]);
      
   if i > 1 && errDNEAq(i) - errDNEAq(i-1) > 0
      disp(['No improvement in estimating dx, therefore exiting. Increase is: ' num2str(errDNEAq(i-1)-errDNEAq(i))]);
      if errDNEAq(i) > 0.02
         disp(['Estimation of dx did not converged! Error was ' num2str(norm([eq; edq])) ' is ' num2str(errDNEAq(i))]);
         res = 1;
         return;
      else
         disp(['Estimation of dx converged. Error was ' num2str(norm([eq; edq])) ' is ' num2str(errDNEAq(i))]);
         return;
      end
   end
   
   q      = randn(NB,1)*0.5;
   dq     = randn(NB,1)*0.5;
   
   dx     = dx + myDNEA.x - x_pri;
   Sx     = myDNEA.Sx + Q;
end

if errdx(end) > 0.02
   disp(['Estimation of dx did not converged! Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(i))]);
   res = 1;
   return;
else
   disp(['Estimation of dx converged. Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(i))]);
   return;
end


end





