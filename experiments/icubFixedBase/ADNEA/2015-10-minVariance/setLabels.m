function setLabels(typeLabel, columns, XY, invertLabels)

old_columns = get(gca, [XY 'Tick']);
old_labels  = get(gca, [XY 'TickLabel']);

[m,n] = size(columns);
if m~=1 && n~=1
   error('Please provide a vector of columns either nx1 or 1xn')
else
   NB = max([m,n]);
end
   
labels = cell(NB, 1);
for i = 1 : NB
   if nargin == 3
      labels{i,1} = [typeLabel num2str(i)];
   else
      labels{i,1} = [typeLabel num2str(NB-i+1)];
   end
end
set(gca, [XY 'Tick'], [old_columns columns], [XY 'TickLabel'], [old_labels; labels])