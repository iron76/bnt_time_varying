function obj = initdBYsubmatrixIndices(obj)

%% Indices for the submatrix
obj.idby = zeros(obj.IDsens.sensorsParams.ny,1);
obj.jdby = zeros(2*obj.IDmodel.n,1);

for i = 1 : obj.IDsens.sensorsParams.ny
   obj.idby(i,1) = obj.IDsens.sensorsParams.sizes{i,1};
end

for i = 1 : obj.IDmodel.n
   obj.jdby((i-1)*2+1:2*i,1) = [1 1]';
end

%% Indeces for the sparse matrix

obj.idby_s = [];
obj.jdby_s = [];

for j = 1 : obj.IDsens.sensorsParams.ny
   sens_str = obj.IDsens.sensorsParams.labels{j,1};
   if length(sens_str) > 7 && strcmp(sens_str(1:7), 'y_omega')
      i = str2double(sens_str(8:end));
      
      for h = 1 : obj.IDstate.n
         if obj.ant(h,i) == 1
            parenth = obj.IDmodel.modelParams.parent(h);
            
            if( parenth ~= 0 )
               obj.idby_s = [obj.idby_s; j];
               obj.jdby_s = [obj.jdby_s; h];
            end
            
            obj.idby_s = [obj.idby_s; j];
            obj.jdby_s = [obj.jdby_s; h+obj.IDstate.n];
         end
      end
   end
end

                 

obj.dby_s = submatrixSparse(obj.idby, obj.jdby, obj.idby_s, obj.jdby_s);