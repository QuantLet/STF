function retval=ultra(x)
% Ultrametric distance between time series.
% x - time series matrix
% 
% Author: Janusz Miśkiewicz, email: jamis@ift.uni.wroc.pl
  [h,k]=size(x);
  retval= sqrt(abs(0.5*(ones(k)-corr(x))));
end