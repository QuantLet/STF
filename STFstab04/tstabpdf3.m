function y=tstabpdf3(x,alpha,sigma,lambda,beta);
x=x(:);%dopisa³em 13-06
if nargin<5
    beta=0;
end;
N=length(x);
if (alpha>0)&&(alpha<=2)&&(alpha~=1)&&(lambda>0)&&(sigma>0)&&(beta<=1)&&(beta>=-1)
    xmax=max(abs(x))+1;
    [xp,yp]=tstabpdf(alpha,sigma,lambda,beta,xmax);
    
%a to przenios³em dla innej interpolacji, ale niepotrzebnie
%     y=zeros(N,1);
%     for i=1:N
%         ind=sum(x(i)>xp);
%         y(i)=yp(ind)+(yp(ind+1)-yp(ind))/(xp(ind+1)-xp(ind))*(x(i)-xp(ind));
%     end;
    
    y=max(interp1q(xp',yp',x),eps);
else
    y=NaN(N,1);
end;