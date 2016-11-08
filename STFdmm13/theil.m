function retval=theil(x,n)
% Converts the given time series into time series of Theil index
% n - size of moving window
% x - analysed time series
%
% Author: Janusz Miśkiewicz, email: jamis@ift.uni.wroc.pl
[n_pocz,k]=size(x);
rozm=n_pocz-n+1;
retval=zeros(rozm,k);
    for i=1:rozm
    	sr=mean(x(i:i+n-1,:));
    	temp=x(i:i+n-1,:)./(ones(n,1)*sr);
    	temp=temp.*log(temp);
    	retval(i,:)=mean(temp);
    end;
end