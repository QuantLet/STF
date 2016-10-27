function [y] = quantilelines(data,step,perc)
if (size(step,1)~=1||step<=0)
  error('quantiles: step must be a positive scalar.');
end
if (size(perc,2)~=1||size(perc,1)==1)
  error('quantiles: perc must be n x 1 vector.');
end

  if(exist('perc')==0)
    perc =(1:9)/10;
  end
  N = size(data,2);
  R = size(data,1);
  begin  = data(1,1,1);
  theend = data(R,1,1);
  numofpoints = (theend-begin)/step+1;
  last        = begin+(numofpoints-1)*step;
  vecstep     = (begin:step:last)';
  y           = zeros(1,size(perc,1));
  i = 1;
  while(i<=numofpoints)
    j = 1;
    vecval = 0;
    while(j<=N)
      aux1 = data(:,j,1);
      aux2 = data(:,j,2);
      pos  = sum(aux1<=vecstep(i));
      if(pos<R)
        vecval=[vecval;aux2(pos)+(vecstep(i)-aux1(pos))*(aux2(pos+1)-aux2(pos))/(aux1(pos+1)-aux1(pos))];
      else
        vecval=[vecval;aux2(pos)];
      end
      j = j + 1;
    end
    y = [y;(quantile(vecval(2:N+1),perc))'];
    i = i + 1;
  end
  y = [vecstep,y(2:numofpoints+1,:)];
end

