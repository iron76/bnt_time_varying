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
      dvdx %% Bidimensional cell array of mdl.n \times 2*mdl.n elements
      % dvdx{i,h} contains \frac{\partial v_i}{\partial q_h} if h <= n
      % otherwise it contains \frac{\partial v_1}{\partial \dot{q}_{h-n}}
      %
      dXupdq %% Cell array of mdl.n elements.  dXupdq{i} contains
      %  \frac{\partial {}^{i}X_{\lambda(i)}}{\partial q_i}
      dDb %% Matrix representing the derivative of D d + b with respect 
      %  to x = [q; dq] \frac{D d + b}{\partial x}
      iDb %% indeces to access the dDb submatrix
      jDb %% indeces to access the dDb submatrix
      
   end
   
   methods
      function b = derivativeIDsolver(mdl,sns)
         if nargin == 0
            error(['You should provide a featherstone-like ' ...
               'model to instantiate stochasticIDsolver'] )
         end
         b = b@stochasticIDsolver(mdl,sns);
         
         %% Initialize own properties
         b.dXupdq  = cell(mdl.n, 1);
         b.dvdx    = cell(mdl.n, 2*mdl.n);
         for i = 1 : mdl.n
            b.dXupdq{i} = zeros(6,6);
            for j = 1 : mdl.n
               b.dvdx{i,j} = zeros(6,1);
               b.dvdx{i,j+mdl.n} = zeros(6,1);
            end
         end
         
         %% Initialize sumatrices
         b = initdDsubmatrixIndices(b);
         b = initdDsubmatrix(b);
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
      end % derivativeIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('derivativeIDsolver disp to be implemented! \n')
      end % disp
      
      function obj = setD(obj,d)
         %% Compute \frac{\partial (Dd+b)}{\partial x} matrix
         [m,n] = size(d);
         if (m ~= obj.IDmodel.modelParams.NB * 26) || (n ~= 1)
            error('[ERROR] The input d should be provided as a column vector with 26*model.NB rows');
         end
         %  as defined in the IJRR Paper
         obj = updateStateDerivativeSubMatrix(obj, d);
      end
      
      function y = simY(obj, d)
         % fprintf('Calling the stochasticIDsolver simY method \n');
         y = cell2mat(obj.IDsens.sensorsParams.Y)*d; % + chol(inv(obj.IDsens.sensorsParams.Sy_inv))*randn(obj.IDmeas.m, 1)
      end
      
   end % methods
end % classdef

