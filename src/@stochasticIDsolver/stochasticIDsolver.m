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
         b   = initSubMatrixIndices(b);
         b   = initSubMatrix(b);

      end % stochasticIDsolver
      
      function disp(b)
         % Display stochasticIDsolver object
         disp@deterministicIDsolver(b)
         fprintf('stochasticIDsolver disp to be implemented! \n')
      end % disp
      
      obj = initSubMatrixIndices(obj);
      obj = initSubMatrix(obj);
   end % methods   
end % classdef

