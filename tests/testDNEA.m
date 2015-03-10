function res = testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel, ymodel, S_dmodel)

% This function checks if the differntial procedure for estimating [d; x]
% from [d_bar; x_bar] works properly. Roughly speaking we start from the
% following equations:
%
% D(x) d + b(x) = 0
% Y(x) d - y    = 0
%
% and the differential estimation should give [d; x] which staisfy the
% above non-linear equations. The estimation is similar to a Newton-like
% optimization step, starting from an initial estimation [d_bar; x_bar].
% The first order approximation of the above non-linear equations around
% [d_bar; x_bar] is given by:
%
% D(x_bar) d + b(x_bar) + dDb(d_bar, x_bar) (x-x_bar) = 0
% Y(x_bar) d - y        + dY(x_bar)         (x-x_bar) = 0
%
% where dDb(d_bar, x_bar) is the derivative of D(x)d+b with respect to x
% evaluated at [d_bar; x_bar] and dY(x_bar) is the derivative of Y(x) with
% respect to x evaluated at x_bar. In order to check the procedure, we
% start by assuming [d_bar; x_bar] such that it satisfies the non-linear
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
% D(x_bar) d + b(x_bar) + dDb(d_bar, x_bar) (x-x_bar) = 0
% Y(x_bar) d - y                                      = 0
%
% The estimation for [d; x] is obtained by solving the above equations in
% [d;x], thus it consits of solving a linear system. If the matrix:
%
% | D(x_bar)    dDb(d_bar, x_bar) |
% | Y(x_bar)    0                 |
%
% is full-rank, a unique solution exists and should coincide with
% [d; x] = [d_bar; x_bar]. When it is not full-rank the function tries to
% modify the sensor placement in a recursive fashion. Heuristically, there
% seems to be quite a good probability of finding a sensor distribution
% which is associated to a full-rank matrix.


ymodel_RNEA = autoSensRNEA(dmodel);
ymodel_RNEA = autoSensStochastic(ymodel_RNEA, 1e-5);

res = 0;
NB  = dmodel.NB;

q         = randn(NB,1)*0.5;
dq        = randn(NB,1)*0.5;
Sx        = diag(1e8*S_dmodel*rand(NB*2, 1));
eq        = 0.05*randn(size(q));
edq       = 0.05*randn(size(dq));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_RNEA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.solveID();
d         = myRNEA.d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_DNEA);

myDNEA = DNEA(myModel, mySens);
myDNEA = myDNEA.setState(q,dq);
y      = myDNEA.simY(d, q, dq);
myDNEA = myDNEA.setState(q+eq,dq+edq);
myDNEA = myDNEA.setY(y);
myDNEA = myDNEA.setStateVariance(Sx);
d_bar  = d + randn(size(d)).*abs(d).*0.1;
myDNEA = myDNEA.setD(d_bar);
myDNEA = myDNEA.solveID();
d1     = myDNEA.d;
dD1    = myDNEA.dDb_s.matrix;
myDNEA = myDNEA.setD(d_bar);
myDNEA = myDNEA.solveID();
d2     = myDNEA.d;
dD2    = myDNEA.dDb_s.matrix;

if (norm(d1-d2) ~= 0)
   res = 1;
   disp('Something wrong with d computed in repetitive calls')
end

if (norm(full(dD1-dD2)) ~= 0)
   res = 1;
   disp('Something wrong with dD computed in repetitive calls')
end


D   = sparse(myDNEA.iDs, myDNEA.jDs, myDNEA.Ds, 19*dmodel.NB, 26*dmodel.NB);
dDb = myDNEA.dDb_s.matrix;
Ys  = myDNEA.IDsens.sensorsParams.Ys;
Ys(:,(end-2*dmodel.NB+1):end) = Ys(:,(end-2*dmodel.NB+1):end) + myDNEA.dby_s.matrix;
DY  = [D dDb; Ys];

if rank(full(DY)) < 26*dmodel.NB + 2*dmodel.NB
   disp([ 'The extended matrix [D;Y] is not full rank! Rank is: ', num2str(rank(full(DY))), ' should be ', num2str(26*dmodel.NB + 2*dmodel.NB)]);
   res = 1;
else
   
   errd = norm(d - d_bar);
   errDNEAd = norm(myDNEA.d - d);
   
   errDNEAq  = norm(myDNEA.x(   1:NB ,1) - q);
   errDNEAdq = norm(myDNEA.x(NB+1:end,1) - dq);
   
   errq  = norm(eq);
   errdq = norm(edq);
   
   disp(['Diff between d.DNEA and d was ' num2str(errd) ' is ' num2str(errDNEAd)]);
   %       if errDNEAd > errd
   %          disp('Result is excessively inaccurate. Test is declared failed!');
   %          res = 1;
   %       end
   
   disp(['Diff between q.DNEA and q was ' num2str(errq) ' is ' num2str(errDNEAq)]);
   if errDNEAq > errq
      disp('Result is excessively inaccurate. Test is declared failed!');
      res = 1;
   end
   
   disp(['Diff between dq.DNEA and dq was ' num2str(errdq) ' is ' num2str(errDNEAdq)]);
   if errDNEAdq > errdq
      disp('Result is excessively inaccurate. Test is declared failed!');
      res = 1;
   end
   
   
   figure
   subplot(211)
   shadedErrorBar(1:26*NB, myDNEA.d, sqrt(diag(myDNEA.Sd)), {'r' , 'LineWidth', 2}, 0);
   hold on
   plot(d)
   
   subplot(212)
   shadedErrorBar(1:2*NB, myDNEA.x, sqrt(diag(myDNEA.Sx)), {'r' , 'LineWidth', 2}, 0);
   hold on
   plot([q; dq])
   
   
end





