function a = initSsubmatrixIndices(a)

a.iSd = zeros(2*a.IDmodel.n,1);
a.jSd = zeros(2*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iSd((i-1)*2+1 : 2*i, 1) = [6 a.IDmodel.jn(i)]';
   a.jSd((i-1)*2+1 : 2*i, 1) = [6 a.IDmodel.jn(i)]';
end
