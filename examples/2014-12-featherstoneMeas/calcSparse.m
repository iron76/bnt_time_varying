function sparseModel = calcSparse(model)

NB  = model.NB;

nb = 0;
nD = zeros(NB, 1);
for i = 1 : NB
  % b1  = Xup{i}*(-a_grav); 
  % OR
  % b1 = crm(v{i})*vJ;
  nb = nb + 6;
  % b2 = crf(v{i})*model.I{i}*v{i};
  nb = nb + 6;
  
  % Dii = [-eye(6) zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) S{i}];
  nD(i) = nD(i) + 6 + jn{i}*6;
  % Dii = [model.I{i} -eye(6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})];
  nD(i) = nD(i) + 36 + 6;
  % Dii = [zeros(6,6) eye(6) -eye(6) zeros(6, jn{i}) -inv(Xa{i}') zeros(6, jn{i})];
  nD(i) = nD(i) + 6 + 6 + 36; 
  % Dii = [zeros(jn{i}, 6) zeros(jn{i}, 6) S{i}' -eye(jn{i}) zeros(jn{i}, 6) zeros(jn{i}, jn{i})];
  nD(i) = nD(i) + 6*jn{i} + jn{i};   
  % Dij = [ Xup{i} zeros(6,6) zeros(6,6) zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i}) zeros(12+jn{i}, 24+2*jn{i})];
  if model.parent(i) ~= 0
    nD(i) = nD(i) + 36;
  end
  
  ind_j  = find(model.parent == i);
  for j = ind_j
    % Dc{i,j} = [ zeros(12, 24+2*jn{i})
    %     zeros(6,6) zeros(6,6) Xup{j}' zeros(6, jn{i}) zeros(6,6) zeros(6, jn{i})
    %     zeros(jn{i}, 24+2*jn{i})];
    nD(i) = nD(i) + 36;
  end
end

pD  = cumsum(nD); 

ib  = zeros(nb, 1);
b1s = zeros(nb, 1);


iD  = zeros(sum(nD), 1);
jD  = zeros(sum(nD), 1);
D1s = zeros(sum(nD), 1);


for i = 1 : NB
  for j = 1 : NB
    inD = [(j-1)*19 (i-1)*19];

    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(1:6));
    aa = aa';   bb = bb';   sparseModel.iD11(1:36,i,j) = aa(:);   sparseModel.jD11(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD12(1:36,i,j) = aa(:);   jD12(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD13(1:36,i,j) = aa(:);   jD13(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+19);
    aa = aa';   bb = bb';   iD14(1: 6,i,j) = aa(:);   jD14(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD21(1:36,i,j) = aa(:);   jD21(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD22(1:36,i,j) = aa(:);   jD22(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD23(1:36,i,j) = aa(:);   jD23(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+19);
    aa = aa';   bb = bb';   iD24(1: 6,i,j) = aa(:);   jD24(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD31(1:36,i,j) = aa(:);   jD31(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(7:12));
    aa = aa';   bb = bb';   iD32(1:36,i,j) = aa(:);   jD32(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(13:18));
    aa = aa';   bb = bb';   iD33(1:36,i,j) = aa(:);   jD33(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+19);
    aa = aa';   bb = bb';   iD34(1: 6,i,j) = aa(:);   jD34(1: 6,i,j) = bb(:);
    
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(1:6));
    aa = aa';   bb = bb';   iD41(1: 6,i,j) = aa(:);   jD41(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(7:12));
    aa = aa';   bb = bb';   iD42(1: 6,i,j) = aa(:);   jD42(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(13:18));
    aa = aa';   bb = bb';   iD43(1: 6,i,j) = aa(:);   jD43(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+19);
    aa = aa';   bb = bb';   iD44(1: 1,i,j) = aa;      jD44(1: 1,i,j) = bb(:);

    inD = [(j-1)*19 19*NB+(i-1)*7];
    
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD15(1:36,i,j) = aa(:);   jD15(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(1:6),inD(2)+7);
    aa = aa';   bb = bb';   iD16(1: 6,i,j) = aa(:);   jD16(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD25(1:36,i,j) = aa(:);   jD25(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(7:12),inD(2)+7);
    aa = aa';   bb = bb';   iD26(1: 6,i,j) = aa(:);   jD26(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+(1:6));
    aa = aa';   bb = bb';   iD35(1:36,i,j) = aa(:);   jD35(1:36,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+(13:18),inD(2)+7);
    aa = aa';   bb = bb';   iD36(1: 6,i,j) = aa(:);   jD36(1: 6,i,j) = bb(:);

    [aa, bb] = meshgrid(inD(1)+19,inD(2)+(1:6));
    aa = aa';   bb = bb';   iD45(1: 6,i,j) = aa(:);   jD45(1: 6,i,j) = bb(:);
    [aa, bb] = meshgrid(inD(1)+19,inD(2)+7);
    aa = aa';   bb = bb';   iD46(1: 1,i,j) = aa(:);   jD46(1: 1,i,j) = bb(:);    
  end
end

pD = [1; pD+1];
