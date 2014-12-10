function [y, Y] = calcMeas(y, a, fB, f, d2q, fx, ny, tau, NB)

Y = cell(ny, NB);
y.values = cell(ny, 1);
for i = 1 : ny
  for j = 1 : NB
    my = y.sizes{i,1};
    Y{i,j} = zeros(my, 26);
    if strcmp(y.labels{i,1}, ['a' num2str(j)])
      Y{i,j} = [eye(6) zeros(6,6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      y.values{i,1} = a{j} + rand(size(a{j}));
    end
    if strcmp(y.labels{i,1}, ['fB' num2str(j)])
      Y{i,j} = [zeros(6,6) eye(6) zeros(6,6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      y.values{i,1} = fB{j}+rand(size(fB{j}));
    end
    if strcmp(y.labels{i,1}, ['f' num2str(j)])
      Y{i,j} = [zeros(6,6) zeros(6,6) eye(6) zeros(6, 1) zeros(6,6) zeros(6, 1)];
      y.values{i,1} = f{j}+rand(size(f{j}));
    end
    if strcmp(y.labels{i,1}, ['tau' num2str(j)])
      Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) eye(1, 1) zeros(1,6) zeros(1, 1)];
      y.values{i,1} = tau(j)+rand(size(tau(j)));
    end
    if strcmp(y.labels{i,1}, ['fx' num2str(j)])
      Y{i,j} = [zeros(6,6) zeros(6,6) zeros(6, 6) zeros(6,1) eye(6) zeros(6, 1)];
      y.values{i,1} = fx{j}+rand(size(fx{j}));
    end
    if strcmp(y.labels{i,1}, ['d2q' num2str(j)])
      Y{i,j} = [zeros(1,6) zeros(1,6) zeros(1,6) zeros(1, 1) zeros(1,6) eye(1, 1)];
      y.values{i,1} = d2q(j)+rand(size(d2q(j)));
    end
    y.values{i,1} = y.values{i,1};
  end
end
