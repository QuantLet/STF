function y=getQnumber(x)
% GETQNUMBER Auxiliary function for STF2loss10
%      y=getQnumber(x) returns the quarterly number of losses. 
         
maxx = max(x);
l = length(x);
y = zeros(maxx,1);
i = 1;
q = 1;

while (i<=l)
  if (x(i)==q)
      y(q) = y(q)+1;
      i = i+1;
  else
      q = q+1;
  end
end
  

