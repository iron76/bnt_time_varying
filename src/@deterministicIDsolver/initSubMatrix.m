function a   = initSubMatrix(a)

a.D = submatrix(a.iD, a.jD);
a.b = submatrix(a.iD, 1);
% Set constant fields for Di,i
for i = 1:a.IDmodel.n
   J = (i-1)*6;
   I = (i-1)*4;
   a.D = set(a.D, -eye(size(a.D(I+1,J+1))),       I+1, J+1);
   a.D = set(a.D, a.IDmodel.S{i},                 I+1, J+6);
   a.D = set(a.D, a.IDmodel.modelParams.I{i},     I+2, J+1);
   a.D = set(a.D, -eye(size(a.D(I+2,J+2))),       I+2, J+2);
   a.D = set(a.D,  eye(size(a.D(I+3,J+2))),       I+3, J+2);
   a.D = set(a.D, -eye(size(a.D(I+3,J+3))),       I+3, J+3);
   a.D = set(a.D, a.IDmodel.S{i}'         ,       I+4, J+3);
   a.D = set(a.D, -eye(size(a.D(I+4,J+4))),       I+4, J+4);
end
