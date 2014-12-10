function [ model ] = autoSensRNEA( n )
%AUTOSENS Generates a random sensor distribution articulated rigid body.
%   Detailed explanation goes here

NB = n;

ny = 0;
for i = 1 : NB
   ny = ny + 1;
   y.sizes{ny,1} = 6;
   y.labels{ny,1} = ['fx'  num2str(i)];
   ny = ny + 1;
   y.sizes{ny,1} = 1;
   y.labels{ny,1} = ['d2q' num2str(i)];
end

model.y = y;
model.m = sum(cell2mat(y.sizes));

Y = cell(ny, NB);
y.values = cell(ny, 1);
for i = 1 : ny
   for j = 1 : NB
      my = y.sizes{i,1};
      Y{i,j} = zeros(my, 26);
      if strcmp(y.labels{i,1}, ['a' num2str(j)])
         Y{i,j} = [eye(6) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      end
      if strcmp(y.labels{i,1}, ['fB' num2str(j)])
         Y{i,j} = [zeros(6,6) eye(6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      end
      if strcmp(y.labels{i,1}, ['f' num2str(j)])
         Y{i,j} = [zeros(6,6) zeros(6,6) eye(6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      end
      if strcmp(y.labels{i,1}, ['tau' num2str(j)])
         Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) eye(1, 1) zeros(1,6) zeros(1, 1)];
      end
      if strcmp(y.labels{i,1}, ['fx' num2str(j)])
         Y{i,j} = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
      end
      if strcmp(y.labels{i,1}, ['d2q' num2str(j)])
         Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) eye(1, 1)];
      end
   end
   model.Y = Y;
end

