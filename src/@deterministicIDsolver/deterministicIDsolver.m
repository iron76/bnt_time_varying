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
      IDmodel, IDsens, d
   end
   
   properties (SetAccess = protected)
      Xup, vJ, v, a, fB, f, Xa, fx, d2q, tau, jn, S
   end
   
   % Class methods
   methods
      function b = deterministicIDsolver(mdl,sns)
         % deterministicIDsolver Constructor function
         if nargin > 0
            if ~checkModel(mdl.modelParams)
               error('You should provide a featherstone-like mdoel')
            end
            b.IDmodel = mdl;
            b.IDsens  = sns;
            b.IDstate = state(mdl.n);
            b.IDmeas  = meas(sns.m);
            b.Xup     = cell(mdl.n, 1);
            b.Xa      = cell(mdl.n, 1);
            b.S       = cell(mdl.n, 1);
            b.vJ      = zeros(6, mdl.n);
            b.d       = zeros(26*mdl.n,1);
            b.v       = zeros(6, mdl.n);
            b.a       = zeros(6, mdl.n);
            b.f       = zeros(6, mdl.n);
            b.fB      = zeros(6, mdl.n);
            b.fx      = zeros(6, mdl.n);
            b.tau     = zeros(mdl.n, 1);
            b.d2q     = zeros(mdl.n, 1);
            b.jn      = zeros(mdl.n, 1);
            for i = 1 : mdl.n
               b.Xup{i}  = zeros(6,6);
               b.Xa{i}   = zeros(6,6);
            end
            for i = 1:mdl.n
               [ ~, b.S{i} ] = jcalc( mdl.modelParams.jtype{i}, 0);
               [~,  b.jn(i)] = size(b.S{i});
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
      
      
      %SETMODELSTATE Set the model position (q) and velocity (dq)
      %   This function sets the position and velocity to be used by the inverse
      %   dynamic solver in following computations.
      %
      %   Genova 6 Dec 2014
      %   Author Francesco Nori
      function obj = setQ(obj,q)
         [n,m] = size(q);
         if (n ~= obj.IDstate.n) || (m ~= 1)
            error('[ERROR] The input q should be provided as a column vector with model.NB rows');
         end
         obj.IDstate.q = q;
         
         for i = 1 : obj.IDstate.n
            [ XJ, ~ ] = jcalc( obj.IDmodel.modelParams.jtype{i}, q(i) );
            obj.Xup{i} = XJ * obj.IDmodel.modelParams.Xtree{i};
         end
      end % Set.q
      
      function obj = setDq(obj,dq)
         [n,m] = size(dq);
         if (n ~= obj.IDstate.n) || (m ~= 1)
            error('[ERROR] The input dq should be provided as a column vector with model.NB rows');
         end
         obj.IDstate.dq = dq;
         
         for i = 1 : obj.IDstate.n
            obj.vJ(:,i) = obj.IDmodel.S(:,i)*dq(i);
         end
         
      end % Set.q
      
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

