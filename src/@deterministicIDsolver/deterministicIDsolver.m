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
      IDmodel, IDsens, d, D, b
   end
   
   properties (SetAccess = protected)
      Xup, vJ, v, a, fB, f, Xa, fx, d2q, tau, jn, iD, jD
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
            a.vJ      = zeros(6, mdl.n);
            a.d       = zeros(26*mdl.n,1);
            a.v       = zeros(6, mdl.n);
            a.a       = zeros(6, mdl.n);
            a.f       = zeros(6, mdl.n);
            a.fB      = zeros(6, mdl.n);
            a.fx      = zeros(6, mdl.n);
            a.tau     = zeros(mdl.n, 1);
            a.d2q     = zeros(mdl.n, 1);
            a.jn      = zeros(mdl.n, 1);
            for i = 1 : mdl.n
               a.Xup{i}  = zeros(6,6);
               a.Xa{i}   = zeros(6,6);
            end
            a.iD = zeros(4*mdl.n,1);
            a.jD = zeros(6*mdl.n,1);
            for i = 1:mdl.n
               a.iD((i-1)*4+1 : 4*i, 1) = [6 6 6 mdl.jn(i)]';
               a.jD((i-1)*6+1 : 6*i, 1) = [6 6 6 mdl.jn(i) 6 mdl.jn(i)]';
               
               % a.hD(i) = 18 +   mdl.jn(i);
               % a.kD(i) = 24 + 2*mdl.jn(i);
            end
            a.D = submatrix(a.iD, a.jD);
            a.b = submatrix(a.iD, 1);
            % Set constant fields for Di,i
            for i = 1:mdl.n
               J = (i-1)*6;
               I = (i-1)*4;
               a.D = set(a.D, -eye(size(a.D(I+1,J+1))), I+1, J+1);
               a.D = set(a.D, mdl.S{i}, I+1,J+6);
               a.D = set(a.D, mdl.modelParams.I{i},     I+2, J+1);
               a.D = set(a.D, -eye(size(a.D(I+2,J+2))), I+2, J+2);
               a.D = set(a.D,  eye(size(a.D(I+3,J+2))), I+3, J+2);
               a.D = set(a.D, -eye(size(a.D(I+3,J+3))), I+3, J+3);
               a.D = set(a.D, mdl.S{i}'               , I+4, J+3);
               a.D = set(a.D, -eye(size(a.D(I+4,J+4))), I+4, J+4);
            end
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
         
         for i = 1 : obj.IDstate.n
            if obj.IDmodel.modelParams.parent(i) ~= 0
               j = obj.IDmodel.modelParams.parent(i);
               I = (i-1)*4;
               J = (j-1)*6;
               obj.D = set(obj.D, obj.Xup{i}, I+1, J+1);
            end
         end
         
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
end % classdef

