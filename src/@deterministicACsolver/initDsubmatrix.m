function a   = initDsubmatrix(a)

a.D = submatrix(a.iD, a.jD);
a.b = submatrix(a.iD, 1);
% Set constant fields for Di,i
for i = 1:a.IDmodel.n
   J = (i-1)*2;
   I = (i-1)*1;
   a.D = set(a.D, -eye(size(a.D(I+1,J+1))),       I+1, J+1);
   a.D = set(a.D, a.IDmodel.S{i},                 I+1, J+2);
end
