function retval=manh(x)
% Manhattan distence between time series normalised by the time series
% length.
% x - time series
%
% Author: Janusz Miśkiewicz, email: jamis@ift.uni.wroc.pl
[h,k]=size(x);
result=zeros(k);
for j=1:k
  for i=1:k
	result(i,j)=abs(mean(x(:,i)-x(:,j)));
  end;
end;
retval= result;
end