function [y, Y, ym, Yx, Yy] = cholMeas(y, ny, NB)

Y = cell(ny, NB);
% y.values = cell(ny, 1);
for i = 1 : ny
  for j = 1 : NB
    my = y.sizes{i,1};
    Y{i,j} = zeros(my, 26);
    if strcmp(y.labels{i,1}, ['a' num2str(j)])
      % Y{i,j} = [eye(6) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      Y{i,j} = sparse(1:6,1:6,ones(6,1), 6, 26);
      y.values{i,1} = sparse(y.values{i,1});
    end
    if strcmp(y.labels{i,1}, ['fB' num2str(j)])
      % Y{i,j} = [zeros(6,6) eye(6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      Y{i,j} = sparse(1:6,7:12,ones(6,1), 6, 26);
      y.values{i,1} = sparse(y.values{i,1});
    end
    if strcmp(y.labels{i,1}, ['f' num2str(j)])
      % Y{i,j} = [zeros(6,6) zeros(6,6) eye(6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      Y{i,j} = sparse(1:6,13:18,ones(6,1), 6, 26);
      y.values{i,1} = sparse(y.values{i,1});
    end
    if strcmp(y.labels{i,1}, ['tau' num2str(j)])
      % Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) eye(1, 1) zeros(1,6) zeros(1, 1)];
      Y{i,j} = sparse(1,19,1,1,26);
      y.values{i,1} = sparse(y.values{i,1});
    end
    if strcmp(y.labels{i,1}, ['fx' num2str(j)])
      % Y{i,j} = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
      Y{i,j} = sparse(1:6,20:25,ones(6,1), 6, 26);
      y.values{i,1} = sparse(y.values{i,1});
    end
    if strcmp(y.labels{i,1}, ['d2q' num2str(j)])
      % Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) eye(1, 1)];
      Y{i,j} = sparse(1,26,1,1,26);
      y.values{i,1} = sparse(y.values{i,1});
    end    
  end
end

yc  = y.values;
ym  = cell2mat(yc);

[ny, ~]   = size(Y);
Yxc = cell(ny,NB);
Yyc = cell(ny,NB);
for i = 1 : ny
  for j = 1 : NB
    Yxc{i,j}  = Y{i,j}(:, 1:19);
    Yyc{i,j}  = Y{i,j}(:, 20:end);
  end
end

Yx = cell2mat(Yxc);
Yy = cell2mat(Yyc);