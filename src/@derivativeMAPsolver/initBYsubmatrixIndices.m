function obj = initBYsubmatrixIndices(obj)

%% Indices for the submatrix
obj.iby = zeros(obj.IDsens.sensorsParams.ny,1);
obj.jby = 1;
for i = 1 : obj.IDsens.sensorsParams.ny
   obj.iby(i,1) = obj.IDsens.sensorsParams.sizes{i,1};
end

%% Indeces for the sparse matrix

obj.iby_s = [];
obj.jby_s = [];

for i = 1 : obj.IDsens.sensorsParams.ny
   sens_str = obj.IDsens.sensorsParams.labels{i,1};
   if length(sens_str)>7 && strcmp(sens_str(1:7), 'y_omega')
      obj.iby_s = [obj.iby_s; i];
      obj.jby_s = [obj.jby_s; 1];
   end
end
obj.by_s = submatrixSparse(obj.iby, obj.jby, obj.iby_s, obj.jby_s);
