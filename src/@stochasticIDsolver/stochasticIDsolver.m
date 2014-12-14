classdef stochasticIDsolver < deterministicIDsolver
   
   properties
      Sd
   end
   
   methods
      function b = stochasticIDsolver(m,y)
         if nargin == 0
            error(['You should provide a featherstone-like ' ... 
             'model to instantiate stochasticIDsolver'] )            
         end
         b = b@deterministicIDsolver(m,y);
      end % stochasticIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('stochasticIDsolver disp to be implemented! \n')
      end % disp
   end % methods
end % classdef

