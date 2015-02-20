function a = initdDsubmatrixIndices(a)

a.iDb = zeros(4*a.IDmodel.n,1);
a.jDb = zeros(2*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iDb((i-1)*4+1 : 4*i, 1) = [6 6 6 a.IDmodel.jn(i)]';
   a.jDb((i-1)*2+1 : 2*i, 1) = [1 1];
   
   % a.hD(i) = 18 +   mdl.jn(i);
   % a.kD(i) = 24 + 2*mdl.jn(i);
end