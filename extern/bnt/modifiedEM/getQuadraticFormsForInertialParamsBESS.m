function [ quadraticForms ] = getQuadraticFormsForInertialParamsBESS( invSigma_fiB , wrench_diag_weights )
%getQuadraticFormsForInertialParamsBESS Get a cell array 10x1 useful to
%compute the E_{t,k}(l(a)^\top \Sigma_{f_i^B}^-1 f_i^B term
quadraticForms = cell(10,1);
for i = 1:10
    quadraticForms{i} = zeros(6,6);
end

for row = 1:6
    for col = 1:6
        acc = zeros(6,1);
        wrench = zeros(6,1);
        acc(row) = 1.0;
        wrench(col) = 1.0;
        buf = (wrench_diag_weights*inertiaRegressor(acc))'*invSigma_fiB*wrench;
        for j = 1:10
            quadraticForms{j}(row,col) = quadraticForms{j}(row,col) + buf(j);
        end
    end
end

end

