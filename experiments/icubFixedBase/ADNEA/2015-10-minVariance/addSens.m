function model = addSens(model, string, S)

% Add a sensor to the sensor model according to the string. Possible
% strings are fi, taui, ai with i = 1 ... model.NB. The associated
% noise variance is in S.

NB = model.NB;
for i = 1 : NB
   if strcmp(['tau' num2str(i)], string)
      s = 1;
      model.ny = model.ny + 1;
      model.labels{end+1} = string;
      model.sizes{end+1} = s;
      model.m = model.m + s;
   end
end

for i = 1 : model.NB
   my = model.sizes{end,1};
   Y{1,i}  = zeros(my, 26);
   Ys{1,i} = zeros(my, 26);
%    if strcmp(model.labels{i,1}, ['a' num2str(j)])
%       model.Y{i,j}  = [eye(6) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%       model.Ys{i,j} = sparse(1:6,1:6,ones(6,1), 6, 26);
%    end
%    if strcmp(model.labels{i,1}, ['fB' num2str(j)])
%       model.Y{i,j}  = [zeros(6,6) eye(6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%       model.Ys{i,j} = sparse(1:6,7:12,ones(6,1), 6, 26);
%    end
%    if strcmp(model.labels{i,1}, ['f' num2str(j)])
%       model.Y{i,j}  = [zeros(6,6) zeros(6,6) eye(6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
%       model.Ys{i,j} = sparse(1:6,13:18,ones(6,1), 6, 26);
%    end
   if strcmp(string, ['tau' num2str(i)])
      Y{1,i}  = [zeros(1,6) zeros(1,6) zeros(1,6) eye(1, 1) zeros(1,6) zeros(1, 1)];
      Ys{1,i} = sparse(1,19,1,1,26);
   end
%    if strcmp(model.labels{i,1}, ['fx' num2str(j)])
%       model.Y{i,j}  = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
%       model.Ys{i,j} = sparse(1:6,20:25,ones(6,1), 6, 26);
%    end
%    if strcmp(model.labels{i,1}, ['d2q' num2str(j)])
%       model.Y{i,j}  = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) eye(1, 1)];
%       model.Ys{i,j} = sparse(1,26,1,1,26);
%    end
end

for i = 1 : model.NB
   Yx{1,i}  = Ys{1,i}(:, 1:19);
   Yy{1,i}  = Ys{1,i}(:, 20:end);
end

model.Ys = [model.Ys; cell2mat([Yx Yy])];
model.Y(end+1, :)  = Y;
      
iSy_s = cell2mat(model.sizes);
jSy_s = cell2mat(model.sizes);
Sy_inv = submatrixSparse(iSy_s, jSy_s, (1:length(iSy_s))', (1:length(jSy_s))');
Sy     = submatrixSparse(iSy_s, jSy_s, (1:length(iSy_s))', (1:length(jSy_s))');

%% Redefine Sy
for i = 1 : model.ny-1
   % model.Sy{i,1} = sMeas.*generateSPDmatrix(dy);
   % S = inv(model.Sy{i,1});
   S = model.Sy(i,i);
   Sy_inv = set(Sy_inv, inv(S), i, i);
   Sy     = set(Sy    ,     S , i, i);
end
Sy_inv = set(Sy_inv, inv(S), model.ny, model.ny);
Sy     = set(Sy    ,     S , model.ny, model.ny);

model.Sy_inv = Sy_inv;
model.Sy     = Sy;
end
      
      
      
      
      