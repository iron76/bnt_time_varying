classdef stochasticIDsolver < deterministicIDsolver
   
   properties
      Sd, Sd_sm, iSd, jSd
   end
   
   methods
      function b = stochasticIDsolver(m,y)
         if nargin == 0
            error(['You should provide a featherstone-like ' ... 
             'model to instantiate stochasticIDsolver'] )            
         end
         b = b@deterministicIDsolver(m,y);
         b   = initVarMatrixIndices(b);
         b   = initVarMatrix(b);

      end % stochasticIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('stochasticIDsolver disp to be implemented! \n')
      end % disp

      function y = simY(obj, d)
         % fprintf('Calling the stochasticIDsolver simY method \n');
         y = cell2mat(obj.IDsens.sensorsParams.Y)*d; % + chol(inv(obj.IDsens.sensorsParams.Sy_inv))*randn(obj.IDmeas.m, 1)
      end 
      
      obj = initVarMatrixIndices(obj);
      obj = initVarMatrix(obj);
   end % methods   
end % classdef

