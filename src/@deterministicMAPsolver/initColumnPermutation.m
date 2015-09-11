function [ a ] = initColumnPermutation( a )
%INITCOLUMNPERMUTATION Summary of this function goes here
%   Detailed explanation goes here
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

end

