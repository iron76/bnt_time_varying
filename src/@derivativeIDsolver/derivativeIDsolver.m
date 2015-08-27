% derivativeIDsolver is a class for wrapping differential inverse dynamics solvers
%
% derivativeIDsolver is a class that is used to wrap mutiple solvers for
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
%     IDmeas - the current measurements: y (class meas)
%    IDmodel - the model of the articulated rigid body (class model)
%     IDsens - the model of the sensor distribution (class sensors)
%          d - dynamic varaibles [a_i, fB_i, f_i, tau_i, fx_i, d2q_i]
%          x - kinematic variables [q, dq]
%         Sx - kinematic vriables variance
%
% METHODS
%   setState - set the current value for x_bar (position q and velocity dq)
%       setY - set the current value for the measurement y
%       setD - set the current value for d_bar
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef derivativeIDsolver < stochasticIDsolver
   properties (SetAccess = protected)
      Sx_inv %% The inverse variance of the current esimation for x
      Sx     %% The variance of the current esimation for x
      dvdx   %% Bidimensional cell array of mdl.n \times 2*mdl.n elements
      % dvdx{i,h} contains \frac{\partial v_i}{\partial q_h} if h <= n
      % otherwise it contains \frac{\partial v_1}{\partial \dot{q}_{h-n}}
      ant    %% ant(h,i) == 1 if h is ancestor of i
      
      dDdq   %% Cell array of mdl.n elements.  dDdq{i} contains
      %  \frac{\partial D(q)}{\partial q_i}
      dXupdq %% Cell array of mdl.n elements.  dXupdq{i} contains
      %  \frac{\partial {}^{i}X_{\lambda(i)}}{\partial q_i}
      dDb %% Matrix representing the derivative of D d + b with respect 
      %  to x = [q; dq] \frac{D d + b}{\partial x}
      iDb %% indeces to access the dDb submatrix
      jDb %% indeces to access the dDb submatrix

      dDb_s %% Matrix representing the derivative of D d + b with respect 
      %  to x = [q; dq] \frac{D d + b}{\partial x}, sparse version
      iDb_s %% indeces to access the dDb submatrix, sparse version
      jDb_s %% indeces to access the dDb submatrix, sparse version  

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
      
      x_bar %% the current estimation for x (around which we linearize)
      d_bar %% the current estimation for d (around which we linearize)
      
      x_prior %% the prior on x
      d_prior %% the prior on d
   end

   properties
      sparsified = 0;     %% check if sparsification was already performed
      S                   %% the sparsification matrix
      x                   %% the result of the estimation for x
   end   
   
   methods
      function b = derivativeIDsolver(mdl,sns)
         if nargin == 0
            error(['You should provide a featherstone-like ' ...
               'model to instantiate stochasticIDsolver'] )
         end
         b = b@stochasticIDsolver(mdl,sns);
         
         %% Initialize own properties
         b.dDdq  = cell(mdl.n, 1);
         b.dXupdq  = cell(mdl.n, 1);
         b.dvdx    = cell(mdl.n, 2*mdl.n);
         for i = 1 : mdl.n
            b.dXupdq{i} = zeros(6,6);
            for j = 1 : mdl.n
               b.dvdx{i,j} = zeros(6,1);
               b.dvdx{i,j+mdl.n} = zeros(6,1);
            end
         end
         
         %% Initialize ancestors
         b.ant = zeros(b.IDstate.n, b.IDstate.n);
         for i = 1 : b.IDstate.n
            h = i;
            while( h ~= 0 )
               parenth = b.IDmodel.modelParams.parent(h);               
               b.ant(h,i) = 1;
               h = parenth;
            end
         end
         
         
         %% Initialize sumatrices
         b = initdDsubmatrixIndices(b);
         b = initdDsubmatrix(b);
         
         b = initBYsubmatrixIndices(b);         
         b = initdBYsubmatrixIndices(b);

      end

      function obj = setState(obj,q,dq)
         obj = setState@stochasticIDsolver(obj,q,dq);
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
         
         %% Value of x around which we linearize        
         obj.x_bar = [q;dq];
      end % derivativeIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('derivativeIDsolver disp to be implemented! \n')
      end % disp

      function obj = setD(obj, d)
         [m,n] = size(d);
         if (m ~= obj.IDmodel.modelParams.NB * 26) || (n ~= 1)
            error('[ERROR] The input d should be provided as a column vector with 26*model.NB rows');
         end
         %  as defined in the IJRR Paper
         obj.d_bar = d;
         obj = updateStateDerivativeSubMatrix(obj, obj.d_bar);
      end
      
      
      function obj = setDprior(obj, d)
         %% Compute \frac{\partial (Dd+b)}{\partial x} matrix
         [m,n] = size(d);
         if (m ~= obj.IDmodel.modelParams.NB * 26) || (n ~= 1)
            error('[ERROR] The input d should be provided as a column vector with 26*model.NB rows');
         end
         %  as defined in the IJRR Paper
         obj.d_prior = d;
      end
      
      function obj = setXprior(obj, x)         
         %% Current prior for the state [q; dq]         
         [m,n] = size(x);
         if (m ~= obj.IDmodel.modelParams.NB * 2) || (n ~= 1)
            error('[ERROR] The input x should be provided as a column vector with 2*model.NB rows');
         end
         %  as defined in the IJRR Paper
         obj.x_prior = x;
      end
      
      function obj = setXvariance(obj,Sx)
         [m,n] = size(Sx);
         if (m ~= obj.IDmodel.modelParams.NB * 2) || (n ~= obj.IDmodel.modelParams.NB * 2)
            error('[ERROR] The input Sx should be a matrix with 2*model.NB rows and columns');
         end
         S = inv(Sx);
         [iS, jS] = find(S);
         obj.Sx_inv = submatrixSparse(ones(m,1), ones(n,1), iS, jS);
         for k = 1 : length(iS)
            obj.Sx_inv = set(obj.Sx_inv, S(iS(k),jS(k)), iS(k), jS(k));
         end
         
         S = Sx;
         [iS, jS] = find(S);
         obj.Sx = submatrixSparse(ones(m,1), ones(n,1), iS, jS);
         for k = 1 : length(iS)
            obj.Sx = set(obj.Sx, S(iS(k),jS(k)), iS(k), jS(k));
         end
      end      
      
      function y = simY(obj, d, q, dq)
         % fprintf('Calling the stochasticIDsolver simY method \n');
         Y = cell2mat(obj.IDsens.sensorsParams.Y);
         y = Y * [d; q; dq] + obj.by_s.matrix; % + chol(inv(obj.IDsens.sensorsParams.Sy_inv))*randn(obj.IDmeas.m, 1)
      end

      function res = eq(obj1, obj2)
         res = 1;
         res = res && isequal(obj1.IDsens.sensorsParams.Sy_inv, obj2.IDsens.sensorsParams.Sy_inv);
         if ~isequal(obj1.IDsens.sensorsParams.Sy_inv, obj2.IDsens.sensorsParams.Sy_inv)
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
         res = res && isequal(obj1.IDmeas.y, obj2.IDmeas.y);
         if ~isequal(obj1.IDmeas.y, obj2.IDmeas.y)
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
      
      %% Compute the derivative of D with respect to qi (i-th element of q)
      function obj = compute_dDdq(obj)
         for i = 1 : obj.IDstate.n
            obj.dDdq{i} = submatrix(obj.iD, obj.jD);
         end
         
         for i = 1 : obj.IDstate.n
            j = obj.IDmodel.modelParams.parent(i);
            if j > 0 
               obj.dDdq{i} = set(obj.dDdq{i}, obj.dXupdq{i}, (i-1)*4+1, (j-1)*6+1);
            end
            
            for j = obj.IDmodel.sparseParams.ind_j{i}
               obj.dDdq{j} = set(obj.dDdq{j}, obj.dXupdq{j}', (i-1)*4+3, (j-1)*6+3);
            end
            
         end
         
      end
      
      
   end % methods
end % classdef

