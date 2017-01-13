function [ ymodel ] = addSens( ymodel , label)
%ADDSENS Add a sensor to a list of sensors

ymodel.ny = ymodel.ny + 1;

for j = 1 : ymodel.NB
   if strcmp(label, ['a' num2str(j)])
      ymodel.m  = ymodel.m + 6;
      ymodel.sizes{ymodel.ny,1} = 6;
      ymodel.labels{ymodel.ny,1} = ['a' num2str(j)];
   end
   if strcmp(label, ['fB' num2str(j)])
      ymodel.m  = ymodel.m + 6;
      ymodel.sizes{ymodel.ny,1} = 6;
      ymodel.labels{ymodel.ny,1} = ['fB' num2str(j)];
   end
   if strcmp(label, ['f' num2str(j)])
      ymodel.m  = ymodel.m + 6;
      ymodel.sizes{ymodel.ny,1} = 6;
      ymodel.labels{ymodel.ny,1} = ['f' num2str(j)];
   end
   if strcmp(label, ['tau' num2str(j)])
      ymodel.m  = ymodel.m + 1;
      ymodel.sizes{ymodel.ny,1} = 1;
      ymodel.labels{ymodel.ny,1} = ['tau' num2str(j)];
   end
   if strcmp(label, ['fx' num2str(j)])
      ymodel.m  = ymodel.m + 6;
      ymodel.sizes{ymodel.ny,1} = 6;
      ymodel.labels{ymodel.ny,1} = ['fx'  num2str(j)];
   end
   if strcmp(label, ['d2q' num2str(j)])
      ymodel.m  = ymodel.m + 1;
      ymodel.sizes{ymodel.ny,1} = 1;
      ymodel.labels{ymodel.ny,1} = ['d2q' num2str(j)];
   end
end

ymodel.Y  = cell(ymodel.ny, ymodel.NB);
ymodel.Ys = cell(ymodel.ny, ymodel.NB);

for i = 1 : ymodel.ny
   for j = 1 : ymodel.NB
      my = ymodel.sizes{i,1};
      ymodel.Y{i,j}  = zeros(my, 26);
      ymodel.Ys{i,j} = zeros(my, 26);
      if strcmp(ymodel.labels{i,1}, ['a' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [diag(d) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,1:6,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['fB' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) diag(d) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,7:12,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['f' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) diag(d) zeros(6, 1) zeros(6,6) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,13:18,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['tau' num2str(j)])
         d = ones(1,1);
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) diag(d) zeros(1,6) zeros(1, 1)];
         ymodel.Ys{i,j} = sparse(1,19,d,1,26);
      end
      if strcmp(ymodel.labels{i,1}, ['fx' num2str(j)])
         d = ones(6,1);
         ymodel.Y{i,j}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) diag(d) zeros(6, 1)];
         ymodel.Ys{i,j} = sparse(1:6,20:25,d, 6, 26);
      end
      if strcmp(ymodel.labels{i,1}, ['d2q' num2str(j)])
         d = ones(1,1);
         ymodel.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) diag(d)];
         ymodel.Ys{i,j} = sparse(1,26,d,1,26);
      end
   end
end

for i = 1 : ymodel.ny
   for j = 1 : ymodel.NB
      Yx{i,j}  = ymodel.Ys{i,j}(:, 1:19);
      Yy{i,j}  = ymodel.Ys{i,j}(:, 20:end);
   end
end
ymodel.Ys = cell2mat([Yx Yy]);
