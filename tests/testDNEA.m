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

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
Sx        = diag(1e8*S_dmodel*rand(dmodel.NB*2, 1));
eq        = 0.0*rand(size(q));
edq       = 0.4*rand(size(dq));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel_RNEA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.solveID();
d_bar     = myRNEA.d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
mySNEA    = mySNEA.setState(q,dq);
y         = mySNEA.simY(d_bar);
mySNEA    = mySNEA.setY(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel_DNEA);
mySens  = sensors(ymodel_DNEA);

myDNEA    = DNEA(myModel, mySens);
myDNEA    = myDNEA.setState(q+eq,dq+edq);
y         = myDNEA.simY([d_bar; q; dq]);
myDNEA    = myDNEA.setY(y);
myDNEA    = myDNEA.setStateVariance(Sx);


if (sum(q-mySNEA.IDstate.q))
   disp('Something wrong with the setQ method');
   res = 1;
end

if (sum(dq-mySNEA.IDstate.dq))
   disp('Something wrong with the setDq method');
   res = 1;
end

if (sum(y-myDNEA.IDmeas.y))
   disp('Something wrong with the setY method');
   res = 1;
end

mySNEA = mySNEA.solveID();

myDNEA = myDNEA.setD(d_bar);
myDNEA = myDNEA.solveID();


D   = sparse(myDNEA.iDs, myDNEA.jDs, myDNEA.Ds, 19*dmodel.NB, 26*dmodel.NB);
dDb = myDNEA.dDb_s.matrix;
Ys  = myDNEA.IDsens.sensorsParams.Ys;
DY  = [D dDb; Ys];

D_SNEA = mySNEA.D.matrix;
b_SNEA = sparse(mySNEA.ibs, ones(size(mySNEA.ibs)), mySNEA.bs, 19*dmodel.NB, 1);
% norm(D_SNEA*mySNEA.d + b_SNEA)

D_DNEA   = myDNEA.D.matrix;
b_DNEA   = sparse(myDNEA.ibs, ones(size(myDNEA.ibs)), myDNEA.bs, 19*dmodel.NB, 1);
% norm(D_DNEA*myDNEA.d + b_DNEA + dDb*(myDNEA.x - myDNEA.x_bar))
% norm(D_DNEA*mySNEA.d + b_DNEA )
% norm(D_DNEA*myDNEA.d + b_DNEA + dDb(:,1:dmodel.NB)*(myDNEA.x(1:dmodel.NB,1) - [q+eq]))


% xx = myDNEA.x;
% norm(D_DNEA*myDNEA.d + b_DNEA - dDb*[q+eq;dq+edq] + dDb(:,1:dmodel.NB)*xx(1:dmodel.NB,1) + dDb(:,dmodel.NB+1:end)*xx(dmodel.NB+1:end,1))
% dqx = pinv(full(dDb(:,dmodel.NB+1:end)))*(-(D_DNEA*myDNEA.d + b_DNEA - dDb*[q+eq;dq+edq] + dDb(:,1:dmodel.NB)*xx(1:dmodel.NB,1)));

if rank(full(DY)) < 26*dmodel.NB + 2*dmodel.NB
   disp([ 'The extended matrix [D;Y] is not full rank! Rank is: ', num2str(rank(full(DY))), ' should be ', num2str(26*dmodel.NB + 2*dmodel.NB)]);
   res = 1;
else
   figure
   subplot(211)
   plot(mySNEA.d)
   hold on
   plot(myDNEA.d, 'r')
   disp(['Diff between d.DNEA and d.SNEA is ' num2str(norm(mySNEA.d-myDNEA.d)/length(mySNEA.d))]);
   if norm(mySNEA.d-myDNEA.d)/length(mySNEA.d) > 0.01*max(mySNEA.d)
      disp('Result is excessively inaccurate. Test is declared failed!');
      res = 1;
   end

   subplot(212)
   plot([q; dq])
   hold on
   plot(myDNEA.x, 'r')
   plot([q+eq;dq+edq], 'g')

   disp(['Diff between x.DNEA and x is ' num2str(norm(myDNEA.x-[q; dq])) ' was ' num2str(norm([eq; edq]))]);
   if norm(myDNEA.x-[q; dq]) > norm([eq; edq])
      disp('Result is excessively inaccurate. Test is declared failed!');
      res = 1;
   end

end





