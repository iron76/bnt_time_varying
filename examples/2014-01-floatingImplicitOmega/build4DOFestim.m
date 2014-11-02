function [h, Sh, Ro, avgTime] = build4DOFestim(index, y, q, dq, o0, do0, p0)


load columnImplicit.mat
defineConst;

c1 = cos(q(3));  c2 = cos(q(4));
s1 = sin(q(3));  s2 = sin(q(4));

% R23
Ro(1:3,1:3,3) = [s1 c1 0;
      -c1  s1 0;
      0    0 1];
% R34
Ro(1:3,1:3,4) = [c2 -s2 0;
      s2  c2 0;
      0    0 1];

Ro(1:3,1:3,5) = [1 0 0;
      0 1 0;
      0 0 1];  
% R01
Ro(1:3,1:3,1) = [1  0 0;
                0  0 1;
                0 -1 0];
% R12  
Ro(1:3,1:3,2) = [0 0 -1;
                0 1  0;
                1 0  0];

j1 = 1:d;
j2 = (1+d*(2-1)):(d*2);
j3 = (1+d*(3-1)):(d*3);
j4 = (1+d*(4-1)):(d*4);
j5 = (1+d*(5-1)):(d*5);

%Useful quantities for implicit NE
param = [l1 l2 m1 m2 r1 r2 g I1z I2z];
xK = [ q'; dq'];
Ci  = CY(xK, param); 
bi  = bY(xK, param);

A = Ci;
b = bi;

indC = 1;
for i = 1 : 4    
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
end
P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
indC = indC+3;
P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
indC = indC+3;
for i = 4: -1 : 1
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
    P(indC:indC+2, indC:indC+2) = sModel.*eye(d);
    indC = indC+3;
end

%fbb
P(indC:indC+2, indC:indC+2) = Sfbb.*eye(d);
b(indC:indC+2) = [0 0 0]';
A(indC:indC+2, colum.fbb) = eye(d);
indC = indC+3;
%ubb
P(indC:indC+2, indC:indC+2) = Subb.*eye(d);
b(indC:indC+2) = [0 0 0]';
A(indC:indC+2, colum.ubb) = eye(d);
indC = indC+3;
%fee
P(indC:indC+2, indC:indC+2) = Sfee.*eye(d);
b(indC:indC+2) = [0 0 0]';
A(indC:indC+2, colum.fee) = eye(d);
indC = indC+3;
%uee
P(indC:indC+2, indC:indC+2) = Suee.*eye(d);
b(indC:indC+2) = [0 0 0]';
A(indC:indC+2, colum.uee) = eye(d);
indC = indC+3;
%d2t
P(indC:indC+3, indC:indC+3) = diag([sUnknown sUnknown sUnknown sUnknown]);
b(indC:indC+3) = [0 0 0 0]';
A(indC:indC+3, colum.d2t) = eye(4);

b  = b  + A(:, [colum.p0 colum.do0])*[p0; do0]; 
A(:, [colum.p0 colum.do0]) = [];

[m, n] = size(P);

tic;
S  = zeros(26, n);
S(1:3  , colum.p(j4)) = eye(3);         S(4:6  , colum.o(j4)) = eye(3);
S(7:9  , colum.f(j2)) = Ro(1:3,1:3,3)'; S(10:12, colum.u(j2)) = Ro(1:3,1:3,3)';
S(13:15, colum.f(j5)) = eye(3);         S(16:18, colum.u(j5)) = eye(3);
S(19   , colum.d2t(3))= 1;              S(20   , colum.d2t(4))= 1;
S(21:23, colum.fbb)   = eye(3);         S(24:26, colum.ubb)   = eye(3);
S(:, [colum.p0 colum.do0]) = [];

R(1:3, 1:3) = Sp4;        R(4:6  , 4:6) = So4;
R(7:9, 7:9) = Sf2;        R(10:12, 10:12) = Su2;
R(13:15, 13:15) = Sf5;        R(16:18, 16:18) = Su5;
R(19,19)= Sd2t3;      R(20   , 20)= Sd2t4;
R(21:23, 21:23)     = Sfbb;       R(24:26, 24:26)   = Subb;

PWinv = inv(A'*inv(P)*A);
Exy = -PWinv*A'*inv(P)*b + PWinv*S'*inv(R+S*PWinv*S')*(y+S*PWinv*A'*inv(P)*b);
Sxy = inv(A'*inv(P)*A + S'*inv(R)*S);
avgTime = toc;

ind = 1;
n = 4;
% [p0(:); do0(:); p(:); do(:); f(:); u(:); d2t(:); c(:); fbb(:); ubb(:); fee(:); uee(:); o(:)]
for i = 1 : n    
    h{index.p(i)}  = Exy(ind:ind+2);
    Sh{index.p(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;
end

for i = 1 : n
    h{index.do(i)}  = Exy(ind:ind+2);
    Sh{index.do(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;
end

for i = 1 : n+1
    h{index.f(i)}  = Exy(ind:ind+2);
    Sh{index.f(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;    
end

for i = 1 : n+1
    h{index.u(i)}  = Exy(ind:ind+2);
    Sh{index.u(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;    
end

for i = 1 : n
    h{index.d2t(i)}  = Exy(ind:ind);
    Sh{index.d2t(i)} = Sxy(ind:ind, ind:ind);
    ind  = ind + 1;    
end

for i = 1 : n
    h{index.c(i)}  = Exy(ind:ind+2);
    Sh{index.c(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;    
end


h{index.fbb}  = Exy(ind:ind+2);
Sh{index.fbb} = Sxy(ind:ind+2, ind:ind+2);
ind  = ind + 3;

h{index.ubb}  = Exy(ind:ind+2);
Sh{index.ubb} = Sxy(ind:ind+2, ind:ind+2);
ind  = ind + 3;

h{index.fee}  = Exy(ind:ind+2);
Sh{index.fee} = Sxy(ind:ind+2, ind:ind+2);
ind  = ind + 3;

h{index.uee}  = Exy(ind:ind+2);
Sh{index.uee} = Sxy(ind:ind+2, ind:ind+2);
ind  = ind + 3;

for i = 1 : n
    h{index.o(i)}  = Exy(ind:ind+2);
    Sh{index.o(i)} = Sxy(ind:ind+2, ind:ind+2);
    ind  = ind + 3;    
end






