classdef stochasticIDsolver < deterministicIDsolver
   
   properties
      S = 0;
   end
   
   methods
      function b = stochasticIDsolver(model)
         if nargin == 0
            error(['You should provide a featherstone-like ' ... 
             'model to instantiate stochasticIDsolver'] )            
         end
         b = b@deterministicIDsolver(model);
      end % stochasticIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('stochasticIDsolver disp to be implemented! \n')
      end % disp
   end % methods
   
   
   methods (Static = true)
      function market_value = calc_value(y)
         market_value = y;
         fprintf('stochasticIDsolver calc_value to be implemented! \n')
      end % calc_value
   end
end % classdef

