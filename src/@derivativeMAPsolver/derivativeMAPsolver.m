% mapMAPsolver is a class for wrapping differential inverse dynamics solvers
%
% derivativeMAPsolver is a class that is used to wrap mutiple solvers for
% computing differential inverse dynamics. Differential inverse dynamic
% solvers compute an estimation of the dynamic (d: forces, torques,
% accelerations, etc.) and kinematic (x: joint position and velocities)
% variables given a set of measurements, possibly redundant. The estimation
% is performed using the following equation:
%
% D(x) d + b(x) = 0         (1)
%
% Y(x) d - y    = 0         (2)
%
% where (1) describes the Newton-Euler eqaution for the articulated rigid
% body and (2) describes the measurement eqaution. The estimation of d and
% x is obtained by first linearizing around [d_bar, x_bar] and then
% computing the maximum a posteriori of the linearized equations. In this
% sense the estimation is nothing but the update step of an implicit Kalman
% filter whose measurement equation is given by (1) and (2). The class
% includes the following properties and methods:
%
% PROPERTIES
%    IDstate - the current articulated rigid body state: q, dq (class state)
%     MAPmeas - the current measurements: y (class meas)
%    IDmodel - the model of the articulated rigid body (class model)
%     MAPsens - the model of the sensor distribution (class sensors)
%          d - dynamic varaibles [a_i, fB_i, f_i, tau_i, fx_i, d2q_i]
%
% METHODS
%   setState - set the current value for x_bar (position q and velocity dq)
%       setY - set the current value for the measurement y
%       setD - set the current value for d_bar
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef derivativeMAPsolver < deterministicMAPsolver
   properties (SetAccess = protected)
      dvdx   %% Bidimensional cell array of mdl.n \times 2*mdl.n elements
      % dvdx{i,h} contains \frac{\partial v_i}{\partial q_h} if h <= n
      % otherwise it contains \frac{\partial v_1}{\partial \dot{q}_{h-n}}
      
      dXupdq %% Cell array of mdl.n elements.  dXupdq{i} contains
      %  \frac{\partial {}^{i}X_{\lambda(i)}}{\partial q_i}

      dbD_s %% Matrix representing the derivative of D d + b with respect 
      %  to x = [q; dq] \frac{D d + b}{\partial x}, sparse version
      ibD_s %% indeces to access the dDb submatrix, sparse version
      jbD_s %% indeces to access the dDb submatrix, sparse version  

      dDd_s %% Matrix representing the derivative of D d + b with respect 
      %  to x = [q; dq] \frac{D d + b}{\partial x}, sparse version
      iDd_s %% indeces to access the dDb submatrix, sparse version
      jDd_s %% indeces to access the dDb submatrix, sparse version  
      
      
      dDdq %% Matrix representing the derivative of D with respect 
      %  to x = q \frac{D}{\partial q}, sparse version. dDdq{i} contains
      %  \frac{\partial D(q)}{\partial q_i}
      iD_s %% indeces to access the dD submatrix, sparse version
      jD_s %% indeces to access the dD submatrix, sparse version        
      
      by_s  %% Matrix representing by(x)
      iby   %% indeces to access the by(x) submatrix
      jby   %% indeces to access the by(x) submatrix
      iby_s %% indeces to access the by(x) submatrix, sparse version
      jby_s %% indeces to access the by(x) submatrix, sparse version  

      dby_s  %% Matrix representing dby(x) derivative of by(x) w.r.t.
      %  x = [q; dq] \frac{by(x)}{\partial x}, sparse version
      idby   %% indeces to access the dby(x) submatrix
      jdby   %% indeces to access the dby(x) submatrix
      idby_s %% indeces to access the dby(x) submatrix, sparse version
      jdby_s %% indeces to access the dby(x) submatrix, sparse version
 end
   
   methods
      function b = derivativeMAPsolver(mdl,sns)
         if nargin == 0
            error(['You should provide a featherstone-like ' ...
               'model to instantiate stochasticMAPsolver'] )
         end
         b = b@deterministicMAPsolver(mdl,sns);
         
         %% Initialize own properties
         b.dDdq  = cell(mdl.n, 1);
         b.dXupdq  = cell(mdl.n, 1);
         b.dvdx    = cell(mdl.n, 2*mdl.n);         
         for i = 1 : mdl.n
            % init dXupdq
            b.dXupdq{i} = zeros(6,6);
            % init dvdx
            for j = 1 : mdl.n
               b.dvdx{i,j} = zeros(6,1);
               b.dvdx{i,j+mdl.n} = zeros(6,1);
            end
            % init dDdq
            b.dDdq{i} = submatrix(b.iD_s, b.jD_s);
         end
         %% Initialize sumatrices         
         b = initBYsubmatrixIndices(b);         
         b = initdBYsubmatrixIndices(b);
         b = initdDsubmatrixIndices(b);
         b = initdDbSubmatrixIndices(b);
      end

      function obj = setState(obj,q,dq)
         obj = setState@deterministicMAPsolver(obj,q,dq);
         %% Compute derivative of the transforms with respect to the parent for all links
         % obj.Xup{i} contains {}^{i}X_{\lambda(i)}
         for i = 1 : obj.IDstate.n
            [ ~, ~, XJderiv ] = jcalcderiv( obj.IDmodel.modelParams.jtype{i}, q(i) );
            obj.dXupdq{i} = XJderiv * obj.IDmodel.modelParams.Xtree{i};
         end
         
         %% Compute derivative of the twist in local frame
         for i = 1 : obj.IDstate.n
            % For each link, compute the derivative against q_h and
            % \dot{q}_h only if and only if h is ancestor of i, all other
            % derivatives are 0
            h = i;
            % X_i_h is {}^iX_h
            X_i_h = eye(6);
            while( h ~= 0 )
               % Position derivative
               parenth = obj.IDmodel.modelParams.parent(h);
               
               % The derivative of the position depend on the velocity of
               % the parent of h, so if the parent of h is the base the
               % derivative is 0
               if( parenth ~= 0 )
                  obj.dvdx{i,h} = X_i_h*obj.dXupdq{h}*obj.v(:,parenth);
               end
               
               % Velocity derivative
               % \frac{\partial v_i}{\partial \dot{q}_h} is equal to
               % {}^iX_h S_h
               obj.dvdx{i,h+obj.IDstate.n} = X_i_h*obj.IDmodel.S{h};
               
               % Now ompute the derivative with respect to the parent
               % Update consequently X_i_h
               h = parenth;
               % Xa{i} = {}^{i}X_{0} => Xa{i}*inv(Xa{h}) = {}^{i}X_{h}
               if  h ~= 0
                  X_i_h = obj.Xa{i}/(obj.Xa{h});
               else
                  X_i_h = obj.Xa{i};
               end
            end
         end
                  
         %% Update the Y and dY matrices         
         obj = updateBYsubmatrix(obj);
         obj = updatedBYsubmatrix(obj);         
      end % derivativeMAPsolver  
      
      
      function y = simY(obj, d, q, dq)
         % fprintf('Calling the stochasticMAPsolver simY method \n');
         Y = cell2mat(obj.MAPsens.sensorsParams.Y);
         y = Y * [d; q; dq] + obj.by_s.matrix; % + chol(inv(obj.MAPsens.sensorsParams.Sy_inv))*randn(obj.MAPmeas.m, 1)
      end

      function res = eq(obj1, obj2)
         res = 1;
         res = res && isequal(obj1.MAPsens.sensorsParams.Sy_inv, obj2.MAPsens.sensorsParams.Sy_inv);
         if ~isequal(obj1.MAPsens.sensorsParams.Sy_inv, obj2.MAPsens.sensorsParams.Sy_inv)
            disp('Sy_inv are different')
         end
         res = res && isequal(obj1.IDmodel.modelParams.Sv_inv.matrix, obj2.IDmodel.modelParams.Sv_inv.matrix);
         if ~isequal(obj1.IDmodel.modelParams.Sv_inv.matrix, obj2.IDmodel.modelParams.Sv_inv.matrix)
            disp('Sv_inv are different')
         end         
         res = res && isequal(obj1.IDmodel.modelParams.Sw_inv.As, obj2.IDmodel.modelParams.Sw_inv.As);
         if ~isequal(obj1.IDmodel.modelParams.Sw_inv.As, obj2.IDmodel.modelParams.Sw_inv.As)
            disp('Sw_inv are different')
         end                  
         res = res && isequal(obj1.Ds, obj2.Ds);
         if ~isequal(obj1.Ds, obj2.Ds)
            disp('Ds are different')
         end
         res = res && isequal(obj1.x_bar, obj2.x_bar);
         if ~isequal(obj1.x_bar, obj2.x_bar)
            disp('x_bar are different')
         end
         res = res && isequal(obj1.d_bar, obj2.d_bar);
         if ~isequal(obj1.d_bar, obj2.d_bar)
            disp('d_bar are different')
         end
         res = res && isequal(obj1.bs, obj2.bs);
         if ~isequal(obj1.bs, obj2.bs)
            disp('bs are different')
         end         
         res = res && isequal(obj1.Sx_inv.As, obj2.Sx_inv.As);
         if ~isequal(obj1.Sx_inv.As, obj2.Sx_inv.As)
            disp('Sx_inv are different')
         end                  
         res = res && isequal(obj1.MAPmeas.y, obj2.MAPmeas.y);
         if ~isequal(obj1.MAPmeas.y, obj2.MAPmeas.y)
            disp('y are different')
         end
         res = res && isequal(obj1.x, obj2.x);
         if ~isequal(obj1.x, obj2.x)
            disp('x are different')
         end
         res = res && isequal(obj1.d, obj2.d);
         if ~isequal(obj1.d, obj2.d)
            disp('d are different')
         end
         
      end
      

      
      %% Compute the derivative of d with respect to q
      function dd_dq = compute_dq(obj, d)
         NB = obj.IDmodel.modelParams.NB;
         
         fwdPerm          = obj.id;
         bckPerm(fwdPerm) = 1:length(fwdPerm);
         
         D = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
         b = sparse(obj.ibs, ones(size(obj.ibs)), obj.bs, 19*NB, 1);
         Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
         Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
         Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;
         Y      = obj.IDsens.sensorsParams.Ys;
         % Y      = Y(1:obj.MAPmeas.m, 1:26*NB);
         
         S_Dinv = Sv_inv;
         S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
         S_Yinv = Sy_inv;
         bY     = zeros(obj.IDmeas.m,1);
         bD     = b;
         
         dbY    = zeros(obj.IDmeas.m, 2*NB);
         
         % Compute dbD
         obj = updatedDbSubmatrix(obj, d);
         dbD = obj.dbD_s.matrix;
         
         % Update dDdq, derivative of D w.r.t. q
         obj = updatedDsubmatrix(obj);
         
         S = inv(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y);
         dd_1  = zeros(26*NB, 2*NB);                     %derivative w.r.t. dq is null
         dd_2  = S*(-Y'*S_Yinv*dbY - D'*S_Dinv*dbD);
         for j = 1 : NB
            dDj   = obj.dDdq{j}.matrix;
            dDj   = dDj(:, bckPerm);
            
            dd_1(:,j) = -S*(D'*S_Dinv*dDj + dDj'*S_Dinv*D)*d(bckPerm,1);
            dd_1(:,j) = dd_1(:,j) + S*(- dDj'*S_Dinv*bD);
         end
         dd_dq = dd_1 + dd_2;
         dd_dq = dd_dq(fwdPerm, :);
         
      end
      
   end % methods
end % classdef

