function [ obj ] = updateDerivativeMatrix( obj )
%UPDATEDERIVATIVEMATRIX Compute \frac{\partial (Dd + b)}{\partial x}
%   Compute the derivative of D(q)d + b(q,\dot{q}) with respect to x
%   (x is defined as (q,\dot{q}).
%   The output is saved in the Ddbx (please suggest a better name) attribute
%   of the input obj.
%   

for h = 1 : obj.IDstate.n

    %% Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b 
    %  with respect to q_h (x_h)
    for i = 1 : obj.IDstate.n
       I = (i-1)*4;
       J = (i-1)*6;

       if i == 1
          obj.b = set(obj.b, obj.Xup{i}*(-obj.IDmodel.g), I+1, 1);
       else
          obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
       end
       obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
       obj.D = set(obj.D, -inv(obj.Xa{i}'), I+3, J+5);
    end

    %% Compute the derivatives of D_{i,j}d_j subvector of Dd + b
    %  with respect to q_h (x_h)
    for i = 1 : obj.IDstate.n
       for j = obj.IDmodel.sparseParams.ind_j{i}
          I = (i-1)*4;
          J = (j-1)*6;
          % f{obj.IDmodel.modelParams.parent(j)} = f{obj.IDmodel.modelParams.parent(j)} + obj.Xup{j}'*f{j};
          % Dc{i,j} = [ zeros(12, 24+2*obj.jn(i))
          %     zeros(6,6) zeros(6,6) obj.Xup{j}' zeros(6, obj.jn(i)) zeros(6,6) zeros(6, obj.jn(i))
          %     zeros(obj.jn(i), 24+2*obj.jn(i))];
          obj.D = set(obj.D, obj.Xup{j}', I+3, J+3);
       end
    end

    %% Compute D_{i,\lambda{i}}d_\lambda{i} subvector of Dd + b 
    %  with respect to q_h (x_h)
    for i = 1 : obj.IDstate.n
       if obj.IDmodel.modelParams.parent(i) ~= 0
          j = obj.IDmodel.modelParams.parent(i);
          I = (i-1)*4;
          J = (j-1)*6;
          obj.D = set(obj.D, obj.Xup{i}, I+1, J+1);
       end
    end
    
    % Only b actually depends on \dot{q}, so we don't have to consider D
    %% Compute the derivatives of D_{i,i}d_i+b_i subvector of Dd + b 
    %  with respect to \dot{q}_h (x_{n+h})
    for i = 1 : obj.IDstate.n
       I = (i-1)*4;
       J = (i-1)*6;

       if i == 1
          obj.b = set(obj.b, obj.Xup{i}*(-obj.IDmodel.g), I+1, 1);
       else
          obj.b = set(obj.b, crm(obj.v(:,i))*obj.vJ(:,i), I+1, 1);
       end
       obj.b = set(obj.b, crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), I+2, 1);
       obj.D = set(obj.D, -inv(obj.Xa{i}'), I+3, J+5);
    end


end

