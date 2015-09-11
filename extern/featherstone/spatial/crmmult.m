function [ vcross ] = crmmult( f, v )
%crfmult computes crf(v)*f


if length(v) == 6

  vcross = [  0    -f(3)  f(2)   0     0     0    ;
	      f(3)  0    -f(1)   0     0     0    ;
	     -f(2)  f(1)  0      0     0     0    ;
	      0    -f(6)  f(5)   0    -f(3)  f(2) ;
	      f(6)  0    -f(4)   f(3)  0    -f(1) ;
	     -f(5)  f(4)  0     -f(2)  f(1)  0 ]*v;


else

  vcross = [  0     0     0    ;
	      f(3)  0    -f(1) ;
	     -f(2)  f(1)  0 ]'*v;
end

end

