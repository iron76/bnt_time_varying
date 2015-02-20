function a = initSsubmatrixIndices(a)

a.iSd = zeros(6*a.IDmodel.n,1);
a.jSd = zeros(6*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iSd((i-1)*6+1 : 6*i, 1) = [6 6 6 a.IDmodel.jn(i) 6 a.IDmodel.jn(i)]';
   a.jSd((i-1)*6+1 : 6*i, 1) = [6 6 6 a.IDmodel.jn(i) 6 a.IDmodel.jn(i)]';
end
