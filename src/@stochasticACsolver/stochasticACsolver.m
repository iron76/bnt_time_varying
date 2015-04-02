classdef stochasticACsolver < deterministicACsolver
   
   properties
      Sd, Sd_sm, iSd, jSd
   end
   
   methods
      function b = stochasticACsolver(m,y)
         if nargin == 0
            error(['You should provide a featherstone-like ' ... 
             'model to instantiate stochasticACsolver'] )            
         end
         b = b@deterministicACsolver(m,y);
         b   = initSsubmatrixIndices(b);
         b   = initSsubmatrix(b);

      end % stochasticACsolver
      
      function disp(b)
         % Display stochasticACsolver object
         disp@deterministicACsolver(b)
         fprintf('stochasticACsolver disp to be implemented! \n')
      end % disp

      function y = simY(obj, d)
         % fprintf('Calling the stochasticACsolver simY method \n');
         y = cell2mat(obj.IDsens.sensorsParams.Y)*d; % + chol(inv(obj.IDsens.sensorsParams.Sy_inv))*randn(obj.IDmeas.m, 1)
      end 
      
      obj = initSsubmatrixIndices(obj);
      obj = initSsubmatrix(obj);
   end % methods   
end % classdef

