function res=samplemef(data, xaxis)
% SAMPLEMEF Sample mean excess function. 
%   RES = SAMPLEMEF (DATA, XAXIS) returns the value of the sample 
%   mean excess function for the vector DATA in the points from vector XAXIS.


sorteddata  = sort (data);
dataLength  =length(sorteddata);
resLength   = length(xaxis);
res         = xaxis;


i=0;
valuesOnTheLeft = 0;
smef = mean (sorteddata(1:dataLength));

while ( i < resLength )
  i = i + 1;
  while (sorteddata(valuesOnTheLeft+1) < xaxis(i))
    valuesOnTheLeft = valuesOnTheLeft + 1;
  end
  smef = mean (sorteddata(valuesOnTheLeft+1:dataLength));
  res(i) = smef - xaxis(i);
end
