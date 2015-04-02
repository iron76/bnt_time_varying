function a = initDsubmatrixIndices(a)

a.iD = zeros(4*a.IDmodel.n,1);
a.jD = zeros(6*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iD((i-1)*4+1 : 4*i, 1) = [6 6 6 a.IDmodel.jn(i)]';
   a.jD((i-1)*6+1 : 6*i, 1) = [6 6 6 a.IDmodel.jn(i) 6 a.IDmodel.jn(i)]';
   
   % a.hD(i) = 18 +   mdl.jn(i);
   % a.kD(i) = 24 + 2*mdl.jn(i);
end

NB = a.IDmodel.n;

id = zeros(6*NB,1);
for i = 1:NB
   id(       (i-1)*4+1 :        4*i, 1) = [6 6 6 a.IDmodel.jn(i)]';
   id(4*NB + (i-1)*2+1 : 4*NB + 2*i, 1) = [6 a.IDmodel.jn(i)]';
end

d_sm = submatrix(id, 1, zeros(26*NB,1));

d_ind = zeros(6*NB,1);
for i = 1 : NB
   d_ind((i-1)*6+1:(i-1)*6+4, 1) =        (i-1)*4+1:       i*4;
   d_ind((i-1)*6+5:(i-1)*6+6, 1) = 4*NB + (i-1)*2+1:4*NB + i*2;
end

[a.id, ~] = indeces(d_sm, d_ind,1);
