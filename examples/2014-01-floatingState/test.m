
param = [l1 l2 m1 m2 r1 r2 g I1z I2z];

x = [0 0 pi/2 -pi/3 100 2 5 6]';

Ci  = CY(x, param);
bi  = bY(x, param);
dC  = dC_dx(x, yi, param);
db  = dC * x;


for i = 1 : 8
    ei =zeros(8,1);
    ei(i) = 1;
    xp = x + ei*1e-12;
    Cip = CY(xp, param);
    bip = bY(xp, param);
    dbi(:,i) = (bip - bi)./1e-12;
    dCi(:,:,i) = (Cip - Ci)./1e-12;
end


ddC = dbi;
for j = 1 : length(yi)
    ddC = ddC + yi(j)*reshape(dCi(:,j,:), 66,8);
end

sum(sum(dbi - db_dx(x, param)))./sum(sum(dbi))
sum(sum(ddC - dC))./sum(sum(ddC))