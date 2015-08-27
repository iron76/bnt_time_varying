function [ quadraticForms ] = getQuadraticFormsForInertialParamsAESS( invSigma_fiB, wrench_weight_diag )
%getQuadraticFormsForInertialParamsBESS Get a cell array 10x10 useful to
%compute the E_{t,k}(l(a)^\top \Sigma_{f_i^B}^-1 l(a) term
quadraticForms = cell(10,10);
for i = 1:10
    for j = 1:10
        quadraticForms{i,j} = zeros(6,6);
    end
end

for row = 1:6
    for col = 1:6
        acc_row = zeros(6,1);
        acc_col = zeros(6,1);
        acc_row(row) = 1.0;
        acc_col(col) = 1.0;
        buf = (wrench_weight_diag*inertiaRegressor(acc_row))'*invSigma_fiB*(wrench_weight_diag*inertiaRegressor(acc_col));
        for i = 1:10
            for j = 1:10
                quadraticForms{i,j}(row,col) = quadraticForms{i,j}(row,col) + buf(i,j);
            end
        end
    end
end

end

