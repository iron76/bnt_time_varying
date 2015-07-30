function [ a_X_b ] = invertAdjointMatrix( b_X_a )
%INVERTADJOINTMATRIX Invert an adjoint matrix
%  Assuming it uses a angular-linear serialization
a_X_b = zeros(6,6);
a_X_b(1:3,1:3) = (b_X_a(1:3,1:3))';
a_X_b(4:6,1:3) = (b_X_a(4:6,1:3))';
a_X_b(4:6,4:6) = (b_X_a(4:6,4:6))';

display(b_X_a*a_X_b)
end

