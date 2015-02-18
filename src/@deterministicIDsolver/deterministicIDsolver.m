% deterministicIDsolver is a class for wrapping inverse dynamics solvers
%
% deterministicIDsolver is a class that is used to wrap mutiple solvers for
% computing inverse dynamics. Inverse dynamic solvers compute an estimation
% of the dynamic varaibles (forces, torques, accelerations, etc.) given a
% set of measurements, possibly redundant. The deterministic solvers, give
% only an estimated value for the dynamic variables, with no variance
% associated to them (see also the class stochasticIDsolver). The class
% include the following properties and methods:
%
% PROPERTIES
%    IDstate - the current articulated rigid body state: q, dq (class state)
%     IDmeas - the current measurements: y (class meas)
%    IDmodel - the model of the articulated rigid body (class model)
%     IDsens - the model of the sensor distribution (class sensors)
%          d - dynamic varaibles [a_i, fB_i, f_i, tau_i, fx_i, d2q_i]
%
% METHODS
%       setQ - set the current value for the position q
%      setDq - set the current value for the velocity dq
%       setY - set the current value for the measurement y
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef deterministicIDsolver
   % file: @deterministicIDsolver/deterministicIDsolver.m
   properties (SetAccess = protected, GetAccess = public)
      IDstate, IDmeas
   end
   
   properties (SetAccess = protected, GetAccess = public)
      IDmodel, IDsens, d, D, b, Ds, bs
      Ddbx %% Matrix of iD rows and 2*n columns, derivative of Dd+q with respect to x
   end
   
   properties (SetAccess = protected)
      Xup %% Cell array of mdl.n elements. Xup{i} contains {}^{i}X_{\lambda(i)}
      Xa  %% Cell array of mdl.n elements. Xup{i} contains {}^{i}X_{0}
      dvdx %% Bidimensional cell array of mdl.n \times 2*mdl.n elements
           % dvdx{i,h} contains \frac{\partial v_i}{\partial q_h} if h <= n
           % otherwise it contains \frac{\partial v_1}{\partial \dot{q}_{h-n}}
           % 
      dXupdq %% Cell array of mdl.n elements.  dXupdq{i} contains 
             %  \frac{\partial {}^{i}X_{\lambda(i)}}{\partial q_i}
      vJ, v, a, fB, f, fx, d2q, tau, iD, jD, iDs, jDs, ibs
   end
   
   % Class methods
   methods
      function a = deterministicIDsolver(mdl,sns)
         % deterministicIDsolver Constructor function
         if nargin > 0
            if ~checkModel(mdl.modelParams)
               error('You should provide a featherstone-like mdoel')
            end
            a.IDmodel = mdl;
            a.IDsens  = sns;
            a.IDstate = state(mdl.n);
            a.IDmeas  = meas(sns.m);
            a.Xup     = cell(mdl.n, 1);
            a.Xa      = cell(mdl.n, 1);
            a.dXupdq  = cell(mdl.n, 1);
            a.dvdx    = cell(mdl.n, 2*mdl.n);
            a.vJ      = zeros(6, mdl.n);
            a.d       = zeros(26*mdl.n,1);
            a.v       = zeros(6, mdl.n);
            a.a       = zeros(6, mdl.n);
            a.f       = zeros(6, mdl.n);
            a.fB      = zeros(6, mdl.n);
            a.fx      = zeros(6, mdl.n);
            a.tau     = zeros(mdl.n, 1);
            a.d2q     = zeros(mdl.n, 1);
            a.Ddbx    = zeros(26*mdl.n, 2*mdl.n);
            for i = 1 : mdl.n
               a.Xup{i}  = zeros(6,6);
               a.Xa{i}   = zeros(6,6);
               a.dXupdq{i} = zeros(6,6);
               for j = 1 : mdl.n
                   a.dvdx{i,j} = zeros(6,1);
                   a.dvdx{i,j+mdl.n} = zeros(6,1);
               end
            end
            a   = initSubMatrixIndices(a);
            a   = initSubMatrix(a);
            
            a   = initSparseMatrixIndices(a);
            a   = initSparseMatrix(a);
            
            % a   = initStateDerivativeSubMatrix(a);
         else
            error(['You should provide a featherstone-like ' ...
               'model to instantiate deterministicIDsolver'] )
         end
      end % deterministicIDsolver
      
      function disp(a)
         % Display a deterministicIDsolver object
         fprintf('deterministicIDsolver disp to be implemented! \n')
         disp(a.IDmodel)
         %fprintf('Description: %s\nDate: %s\nType: %s\nCurrent Value: $%4.2f\n',...
         %   a.Description,a.Date,a.Type,a.CurrentValue);
      end % disp
      
      
      %SETSTATE Set the model position (q) and velocity (dq)
      %   This function sets the position and velocity to be used by the inverse
      %   dynamic solver in following computations.
      %
      %   Genova 6 Dec 2014
      %   Author Francesco Nori
      
      function obj = setState(obj,q,dq)
         [n,m] = size(q);
         if (n ~= obj.IDstate.n) || (m ~= 1)
            error('[ERROR] The input q should be provided as a column vector with model.NB rows');
         end
         obj.IDstate.q = q;
         
         %% Compute transforms with respect to the parent for all links
         % obj.Xup{i} contains {}^{i}X_{\lambda(i)}
         for i = 1 : obj.IDstate.n
            [ XJ, ~ ] = jcalc( obj.IDmodel.modelParams.jtype{i}, q(i) );
            obj.Xup{i} = XJ * obj.IDmodel.modelParams.Xtree{i};
         end
         
         %% Compute derivative of the transforms with respect to the parent for all links
         % obj.Xup{i} contains {}^{i}X_{\lambda(i)}
         for i = 1 : obj.IDstate.n
            [ ~, ~, XJderiv ] = jcalcderiv( obj.IDmodel.modelParams.jtype{i}, q(i) );
            obj.dXupdq{i} = XJderiv * obj.IDmodel.modelParams.Xtree{i};
         end
         
         %% Compute transforms with respect to the base for all links
         % obj.Xa{i} contains {}^{i}X_{0}
         for i = 1:length(obj.IDmodel.modelParams.parent)
            if obj.IDmodel.modelParams.parent(i) == 0
               obj.Xa{i} = obj.Xup{i};
            else
               obj.Xa{i} = obj.Xup{i} * obj.Xa{obj.IDmodel.modelParams.parent(i)};
            end
         end
         
         [n,m] = size(dq);
         if (n ~= obj.IDstate.n) || (m ~= 1)
            error('[ERROR] The input dq should be provided as a column vector with model.NB rows');
         end
         obj.IDstate.dq = dq;
         
         %% Compute twist in local frame
         % obj.v(:,i) contains {}^{i}v_i
         for i = 1 : obj.IDstate.n
            obj.vJ(:,i) = obj.IDmodel.S{i}*dq(i);
            if obj.IDmodel.modelParams.parent(i) == 0
               obj.v(:,i) = obj.vJ(:,i);
            else
               obj.v(:,i) = obj.Xup{i}*obj.v(:,obj.IDmodel.modelParams.parent(i)) + obj.vJ(:,i);
            end
         end
         
         %% Compute derivative of the twist in local frame
         for i = 1 : obj.IDstate.n
             % For each link, compute the derivative against q_h and
             % \dot{q}_h only if and only if h is ancestor of i, all other
             % derivatives are 0
             h = i;
             % X_i_h is {}^iX_h
             X_i_h = eye(6,6);
             while( h ~= 0 )
                 % Position derivative
                 parenth = obj.IDmodel.modelParams.parent(h);
                 
                 % The derivative of the position depend on the velocity of
                 % the parent of h, so if the parent of h is the base the
                 % derivative is 0
                 if( parenth ~= 0 )
                     obj.dvdx{i,h+obj.IDstate.n} = X_i_h*obj.dXupdq{h}*obj.v(:,parenth);
                 end
                 
                 % Velocity derivative
                 % \frac{\partial v_i}{\partial \dot{q}_h} is equal to
                 % {}^iX_h S_h
                 obj.dvdx{i,h+obj.IDstate.n} = X_i_h*obj.IDmodel.S{i};

                 % Now ompute the derivative with respect to the parent
                 % Update consequently X_i_h
                 h = parenth;
                 X_i_h = X_i_h*obj.Xup{i};
             end
         end
         
         %% Compute D matrix and b vector, as defined in the IJRR Paper
         % 
         obj = updateSubMatrix(obj);
         
         %% Compute \frac{\partial (Dd+b)}{\partial x} matrix
         %  as defined in the IJRR Paper
         % obj = updateStateDerivativeSubMatrix(obj);
         
         %% Update the sparse representation of the matrix D
         obj = updateSparseMatrix(obj);

      end % setState
      
      function obj = setY(obj,y)
         [m,n] = size(y);
         if (m ~= obj.IDmeas.m) || (n ~= 1)
            error('[ERROR] The input y should be provided as a column vector with model.NB rows');
         end
         obj.IDmeas.y = y;
      end 
      
      function y = simY(obj, d)
         % fprintf('Calling the deterministicIDsolver simY method \n');
         y = cell2mat(obj.IDsens.sensorsParams.Y)*d;
      end 
      
      obj = solveID(obj)
   end
   
   methods
      obj = initSparseMatrixIndices(obj);
      obj = initSubMatrixIndices(obj);
      obj = initSubMatrix(obj);
      obj = initSparseMatrix(obj);
      % obj = initStateDerivativeSubMatrix(obj);

      
      obj = updateSubMatrix(obj);
      obj = updateSparseMatrix(obj);
      % obj = updateStateDerivativeSubMatrix(obj);
      
   end
end % classdef

