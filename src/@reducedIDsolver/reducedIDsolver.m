% reducedIDsolver is a class for wrapping inverse dynamics solvers
%
% reducedIDsolver is a class that is used to wrap mutiple solvers for
% computing inverse dynamics. Inverse dynamic solvers compute an estimation
% of the dynamic varaibles (forces, torques, accelerations, etc.) given a
% set of measurements, possibly redundant. The reduced solvers, give
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

classdef reducedIDsolver
   % file: @reducedIDsolver/reducedIDsolver.m
   properties (SetAccess = protected, GetAccess = public)
      IDstate, IDmeas
   end
   
   properties (SetAccess = protected, GetAccess = public)
      IDmodel, IDsens, d, D, b, Ds, bs
   end
   
   properties (SetAccess = protected)
      Xup, vJ, v, a, f, Xa, d2q, iD, jD, iDs, jDs, ibs
   end
   
   % Class methods
   methods
      function a = reducedIDsolver(mdl,sns)
         % reducedIDsolver Constructor function
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
            a.vJ      = zeros(6, mdl.n);
            a.d       = zeros(26*mdl.n,1);
            a.v       = zeros(6, mdl.n);
            a.a       = zeros(6, mdl.n);
            a.f       = zeros(6, mdl.n);
            a.d2q     = zeros(mdl.n, 1);
            for i = 1 : mdl.n
               a.Xup{i}  = zeros(6,6);
               a.Xa{i}   = zeros(6,6);
            end
            a   = initSubMatrixIndices(a);
            a   = initSubMatrix(a);
            
            % a   = initSparseMatrixIndices(a);
            % a   = initSparseMatrix(a);               
         else
            error(['You should provide a featherstone-like ' ...
               'model to instantiate reducedIDsolver'] )
         end
      end % reducedIDsolver
      
      function disp(a)
         % Display a reducedIDsolver object
         fprintf('reducedIDsolver disp to be implemented! \n')
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
         
         for i = 1 : obj.IDstate.n
            [ XJ, ~ ] = jcalc( obj.IDmodel.modelParams.jtype{i}, q(i) );
            obj.Xup{i} = XJ * obj.IDmodel.modelParams.Xtree{i};
         end
         
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
         
         for i = 1 : obj.IDstate.n
            obj.vJ(:,i) = obj.IDmodel.S{i}*dq(i);
            if obj.IDmodel.modelParams.parent(i) == 0
               obj.v(:,i) = obj.vJ(:,i);
            else
               obj.v(:,i) = obj.Xup{i}*obj.v(:,obj.IDmodel.modelParams.parent(i)) + obj.vJ(:,i);
            end
         end
         
         obj = updateSubMatrix(obj);
         % obj = updateSparseMatrix(obj);

      end % setState
      
      function obj = setY(obj,y)
         [m,n] = size(y);
         if (m ~= obj.IDmeas.m) || (n ~= 1)
            error('[ERROR] The input y should be provided as a column vector with model.NB rows');
         end
         obj.IDmeas.y = y;
      end % Set.q
      
      obj = solveID(obj)
   end
   
   methods
      % obj = initSparseMatrixIndices(obj);
      % obj = initSparseMatrix(obj);
      % obj = updateSparseMatrix(obj);

      obj = initSubMatrixIndices(obj);
      obj = initSubMatrix(obj);
      obj = updateSubMatrix(obj);
   end
end % classdef

