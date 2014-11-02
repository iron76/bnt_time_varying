function W = commuteOnOrder2(A, B, a, b)
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
[y,ind] = sort([a, b]);
W = [Wt{ind(1)} Wt{ind(2)}];
