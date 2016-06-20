function [P, Q, WL, WR] = factorize(obj)

% Author: Francesco Nori
% Genova, June 2016
%
% Given an forward-dynamics problem formulated as [D; Y]*d + [bD; -y] = 0
% with d organized as follows:
%
%   d   = [dx_1, dx_2, ..., dx_obj.IDstate.n, dy_1, dy_2, ..., dy_obj.IDstate.n]
%
%   where:
%
%   dx_i = [a_i, fB_i, f_i, tau_i]
%
%   dy_i = [fx_i, d2q_i]
%
% computes permutation matrices P and Q such that P*[D;Y]*Q = L with L
% lower triangular. 

NB = obj.IDmodel.modelParams.NB;


for i = 1 : length(obj.jd2q)
   jd2q_a(1+7*(i-1):7*i,1) = [obj.jd2q(i); obj.ja(1+6*(i-1):6*i)];
end
for i = 1 : length(obj.jtau)
   itau_a(1+7*(i-1):7*i,1) = [obj.itau(i); obj.ia(1+6*(i-1):6*i)];
end


p  = [obj.ifB(end:-1:1, 1); obj.iF(end:-1:1, 1); itau_a];
q  = [obj.jfx; obj.jtau; obj.jfB(end:-1:1, 1); obj.jF(end:-1:1, 1); jd2q_a];
D  = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
Y  = obj.IDsens.sensorsParams.Ys;

% Let's build the matrix Ws such that (if lambda_i ~=0):
%
%   pi_a = f_i - Ii_a a_lambda_i
%
% and (if lambda_i == 0)
%
%   pi_A = f_i - Ii_A a_i
%
% In the case lambda_i = i-1, lambda_1 = 0 we have:
%
%   | p1_A |   | f1 |     | a1 |  
%   | p2_A |   | f2 |     | a2 |  
%   |  ..  | = | .. | + W | .. |
%   | pN_A |   | fN |     | aN | 
%
%              | -I1_A     0   ...    0    0  |
%              |     0 -I2_A   ...    0    0  |
%        Ws  = |     0     0   ...    0    0  |
%              |                              |
%              |     0     0   ...    0 -In_A |
%
%
Ws = zeros(6*NB, 6*NB);
Ws_inv = zeros(6*NB, 6*NB);
for i = 1 : NB
   %if obj.IDmodel.modelParams.parent(i) == 0
      Ws(1+6*(i-1):6*i, 1+6*(i-1):6*i) = -obj.IA{i};
      Ws_inv(1+6*(i-1):6*i, 1+6*(i-1):6*i) = obj.IA{i};
   %else
   %   j = obj.IDmodel.modelParams.parent(i);
   %   Ws(1+6*(i-1):6*i, 1+6*(j-1):6*j) = -obj.Ia{i}*obj.Xup{i};
   %end
end

% Let's build the matrix Wr such that:
%
%   fi_b = fi_B - Ii ai
%
%   | f1_b |   | f1_B |      | a1 |  
%   | f2_b |   | f2_B |      | a2 |  
%   |  ..  | = |  ..  | + Wr | .. |
%   | fN_b |   | fN_B |      | aN | 
%
%              | -I1    0   ...   0 |
%              |   0  -I2   ...   0 |
%        Wr  = |   0    0   ...   0 |
%              |                    |
%              |   0    0   ... -IN |
%
%
Wr = zeros(6*NB, 6*NB);
Wr_inv = zeros(6*NB, 6*NB);
for i = 1 : NB
   Wr(1+6*(i-1):6*i, 1+6*(i-1):6*i) = -obj.IDmodel.modelParams.I{i};
   Wr_inv(1+6*(i-1):6*i, 1+6*(i-1):6*i) = obj.IDmodel.modelParams.I{i};
end

W = eye(26*NB, 26*NB);
W(obj.jF, obj.ja) = Ws;
W(obj.jfB, obj.ja) = Wr;
W = sparse(W);

W_inv = eye(26*NB, 26*NB);
W_inv(obj.jF, obj.ja) = Ws_inv;
W_inv(obj.jfB, obj.ja) = Wr_inv;
WR = sparse(W_inv);

WL = eye(19*NB, 19*NB);
for i = 1:NB 
      
   I = obj.ia(1+6*(i-1):6*i);
   J = obj.ja(1+6*(i-1):6*i);
   K = obj.itau(i);
   % D(K,J) = obj.IDmodel.S{i}' * obj.IA{i}
   W1 = zeros(19*NB, 19*NB);
   W1(K,I) = obj.IDmodel.S{i}' * obj.IA{i};
   W1 = eye(19*NB, 19*NB) + W1;
   W1 = sparse(W1);
   WL = W1*WL;
   % D(K,J)

   J = obj.jd2q(i);
   K = obj.itau(i);
   % D(K,J) = obj.IDmodel.S{i}' * obj.IA{i} * obj.IDmodel.S{i}
   W2 = eye(19*NB, 19*NB);
   W2(K, K) = inv(obj.IDmodel.S{i}' * obj.IA{i} * obj.IDmodel.S{i});
   W2 = sparse(W2);
   WL = W2*WL;
   
   I = obj.itau(i);
   J = obj.jd2q(i);
   K = obj.ia(1+6*(i-1):6*i);
   % D(K,J) = obj.IDmodel.S{i}
   W3 = zeros(19*NB, 19*NB);
   W3(K,I) = -obj.IDmodel.S{i};
   W3 = eye(19*NB, 19*NB) + W3;
   W3 = sparse(W3);
   WL = W3*WL;
   %D(obj.itau(2), obj.jd2q(2))
end

for i = 1:NB 
   for j = obj.IDmodel.sparseParams.ind_j{i}
      I = obj.ia(1+6*(j-1):6*j);
      J = obj.ja(1+6*(j-1):6*j);
      K = obj.iF(1+6*(i-1):6*i);
      % D(K,J) = obj.Xup{j}' * obj.IA{j}
      W4 = zeros(19*NB, 19*NB);
      W4(K,I) = obj.Xup{j}' * obj.IA{j};
      W4 = eye(19*NB, 19*NB) + W4;
      W4 = sparse(W4);
      WL = W4*WL;
   end
end

I  = eye(size([D;Y]));
Q  = I(:,q);

ID = eye(size(D*D'));
IY = eye(size(Y*Y'));
P  = blkdiag(IY, ID(p,:))*[0*Y*D', IY; ID, 0*D*Y'];

P = sparse(P);
Q = sparse(Q);

% The resulting P*[WL*D*WR; Y]*Q is lower triangular
return




