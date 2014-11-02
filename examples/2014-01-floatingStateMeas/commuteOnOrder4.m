function W = commuteOnOrder4(A, B, C, D, a, b, c, d)
% if (a < b && b < c)
%     W = [A B C];
% elseif (b < a && a < c)
%     W = [B A C];
% elseif (a < c && c < b)
%     W = [B C A];
% elseif (b < c && c < a)
%     W = [A C B];
% elseif (c < a && a < b)
%     W = [C B A];
% elseif (c < b && b < a)
%     W = [C A B];
% end

Wt{1} = A;
Wt{2} = B;
Wt{3} = C;
Wt{4} = D;
[y,ind] = sort([a, b, c, d]);
W = [Wt{ind(1)} Wt{ind(2)} Wt{ind(3)} Wt{ind(4)}];
