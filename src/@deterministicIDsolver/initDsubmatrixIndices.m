function a = initDsubmatrixIndices(a)

a.iD = zeros(4*a.IDmodel.n,1);
a.jD = zeros(6*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iD((i-1)*4+1 : 4*i, 1) = [6 6 6 a.IDmodel.jn(i)]';
   a.jD((i-1)*6+1 : 6*i, 1) = [6 6 6 a.IDmodel.jn(i) 6 a.IDmodel.jn(i)]';
   
   % a.hD(i) = 18 +   mdl.jn(i);
   % a.kD(i) = 24 + 2*mdl.jn(i);
end
