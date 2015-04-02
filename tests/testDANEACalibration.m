function res = testDANEACalibration(dmodel_DANEA, ymodel_DANEA, dmodel, ymodel, S_dmodel)

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
eq        = 15*randn(NB,1)*pi/180;
edq       = 15*randn(NB,1)*pi/180;
dx        = zeros(2*NB,1);
q         = randn(NB,1)*0.5;
dq        = randn(NB,1)*0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DANEA);
mySens  = sensors(ymodel_RNEA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDANEA       = DANEA(model(dmodel_DANEA), sensors(ymodel_DANEA));
% myModel = myModel.updateVarianceDNEA(0, myDANEA.dDb_s);
% mySens  = mySens.updateVarianceDNEA(0, myDANEA.dby_s);
% myDANEA  = DNEA(myModel, mySens);


for i = 1 : T
   myRNEA    = myRNEA.setState(q,dq);
   myRNEA    = myRNEA.solveID();
   d         = myRNEA.d;
   d_red     = zeros(7*NB,1);
   for j = NB : -1 : 1
      d_red(1+(j-1)*7 : 6+(j-1)*7, 1 ) = d( 1+(j-1)*26 :  6+(j-1)*26, 1 );
      d_red(7+(j-1)*7 : 7+(j-1)*7, 1 ) = d(26+(j-1)*26 : 26+(j-1)*26, 1 );
   end
   % disp(['RNEA Error in D d + b is: ' num2str(norm(myRNEA.D.matrix*d + myRNEA.b.matrix))]);
   
   d_pri  = d_red + randn(size(d_red)).*abs(d_red).*0.1;
   d_bar  = d_pri;
   
   x      = [q; dq];
   x_bar  = [q; dq] + [eq; edq];
   x_pri  = x_bar + dx;
   
   disp(['Distance from real x is: ' num2str(norm(x_pri - x))]);
      
   myDANEA = myDANEA.setState(x(1:NB ,1),x(NB+1:end,1));
   y       = myDANEA.simY(d_red,  x(1:NB ,1),x(NB+1:end,1));
   
   myDANEA = myDANEA.setState(x_pri(1:NB ,1),x_pri(NB+1:end,1));
   myDANEA = myDANEA.setY(y);
   
   myDANEA = myDANEA.setD(d_bar);
   myDANEA = myDANEA.setDprior(d_pri);
   myDANEA = myDANEA.setXprior(x_pri);
   myDANEA = myDANEA.setXvariance(Sx);
   myDANEA = myDANEA.solveID();
      
   errDNEAd(i)  = norm(myDANEA.d - d_pri);
   errDNEAq(i)  = norm(myDANEA.x(   1:NB ,1)  -  q);
   errDNEAdq(i) = norm(myDANEA.x(NB+1:end,1)  - dq);
   
   errd(i)  = norm(d_pri-d_bar);
   errq(i)  = norm(eq);
   errdq(i) = norm(edq);
   errdx(i) = norm(dx+[eq; edq]);
   
   D   = myDANEA.D(1:NB, 1:2*NB);
   dDb = myDANEA.dDb_s.matrix;
   Ys  = myDANEA.IDsens.sensorsParams.Ys;
   Ys(:,(end-2*dmodel.NB+1):end) = Ys(:,(end-2*dmodel.NB+1):end) + myDANEA.dby_s.matrix;
   DY  = [D dDb; Ys];
   %norm(DY*[d_red; x*0]-[zeros(6*NB,1); y])
      
   
   
   if rank(full(DY)) < 7*dmodel.NB + 2*dmodel.NB
      disp([ 'The extended matrix [D;Y] is not full rank! Rank is: ', num2str(rank(full(DY))), ' should be ', num2str(7*dmodel.NB + 2*dmodel.NB)]);
      res = 1;
   end
      
   if i > 1 && errDNEAq(i) - errDNEAq(i-1) > 1
      disp(['No improvement in estimating dx, therefore exiting. Increase is: ' num2str(errDNEAq(i-1)-errDNEAq(i))]);
      if errdx(i) > 0.02
         disp(['Estimation of dx did not converged! Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(i))]);
         res = 1;
         return;
      else
         disp(['Estimation of dx converged. Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(i))]);
         return;
      end
   end

   q      = randn(NB,1)*0.5;
   dq     = randn(NB,1)*0.5;

   dx     = dx + myDANEA.x - x_pri;
   Sx     = myDANEA.Sx + Q;
end

if errdx(end) > 0.02
   disp(['Estimation of dx did not converged! Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(end))]);
   res = 1;
   return;
else
   disp(['Estimation of dx converged. Error was ' num2str(norm([eq; edq])) ' is ' num2str(errdx(end))]);
   return;
end


end
