function a = initDsubmatrixIndices(a)

a.iD = zeros(1*a.IDmodel.n,1);
a.jD = zeros(2*a.IDmodel.n,1);
for i = 1:a.IDmodel.n
   a.iD((i-1)*1+1 : 1*i, 1) = 6;
   a.jD((i-1)*2+1 : 2*i, 1) = [6 a.IDmodel.jn(i)]';
   
   % a.hD(i) = 18 +   mdl.jn(i);
   % a.kD(i) = 24 + 2*mdl.jn(i);
end

NB = a.IDmodel.n;

id = zeros(2*NB,1);
for i = 1:NB
   id(     (i-1)*1+1 :      1*i, 1) = 6;
   id(NB + (i-1)*1+1 : NB + 1*i, 1) = a.IDmodel.jn(i);
end

d_sm = submatrix(id, 1, zeros(7*NB,1));

d_ind = zeros(2*NB,1);
for i = 1 : NB
   d_ind((i-1)*2+1:(i-1)*2+1, 1) =      (i-1)*1+1:     i*1;
   d_ind((i-1)*2+2:(i-1)*2+2, 1) = NB + (i-1)*1+1:NB + i*1;
end

[a.id, ~] = indeces(d_sm, d_ind,1);
