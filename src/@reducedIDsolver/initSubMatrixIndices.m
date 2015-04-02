function a = initSubMatrixIndices(a)

a.iD = zeros(2*a.IDmodel.n,1);
a.jD = zeros(3*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iD((i-1)*4+1 : 4*i, 1) = [6 6]';
   a.jD((i-1)*6+1 : 6*i, 1) = [6 6 a.IDmodel.jn(i)]';
   
   % a.hD(i) = 18 +   mdl.jn(i);
   % a.kD(i) = 24 + 2*mdl.jn(i);
end
