function res = testDNEA(dmodel, ymodel)

res = 0;

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);
Sx        = diag(rand(dmodel.NB*2, 1)*1e-12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mySNEA    = SNEA(myModel, mySens);
mySNEA    = mySNEA.setState(q,dq);
mySNEA    = mySNEA.setY(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myDNEA    = DNEA(myModel, mySens);
myDNEA    = myDNEA.setState(q,dq);
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

if (sum(y-mySNEA.IDmeas.y))
   disp('Something wrong with the setY method');
   res = 1;
end

mySNEA = mySNEA.solveID();

myDNEA = myDNEA.setD(mySNEA.d);
myDNEA = myDNEA.solveID();

D  = sparse(myDNEA.iDs, myDNEA.jDs, myDNEA.Ds, 19*dmodel.NB, 26*dmodel.NB); 
DY = [D myDNEA.dDb_s.matrix; myDNEA.IDsens.sensorsParams.Ys zeros(myDNEA.IDmeas.m, dmodel.NB*2)];
if rank(full(DY)) < 26*dmodel.NB + 2*dmodel.NB
   disp('The extended matrix [D;Y] is not full rank!');
   res = 1;
end   
   
   
disp(['Diff between DNEA and SNEA is ' num2str(norm(mySNEA.d-myDNEA.d))]);
if norm(mySNEA.d-myDNEA.d) > 1
   disp('Result is excessively inaccurate. Test is declared failed!');
   res = 1;
end




