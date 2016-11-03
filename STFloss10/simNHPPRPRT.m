function [t,y]=simNHPPRPRT(u,theta,lambda,parlambda,distrib,params,time,val,T)
% SIMNHPPRPRT Real-life trajectory of the risk process.
%       [t,y]=simNHPPRPRT(u,theta,lambda,parlambda,distrib,params,time,val,
%       T) returns the real-life trajectory of the risk process 
%       with claim sizes specified by VAL and the moments of losses by
%       TIME. The premium corresponds to the non-homogeneous Poisson
%       process with inensity specified by LAMBDA (0 - sine function, 
%       1 - linear function, 3 - sine square function) with parameters in
%       PARLAMBDA and mean claim size from a distribution specified by
%       DISTRIB with parameters in PARAMS. THETA is the relative safety
%       loading, U the initial capital and T the time horizon.

  lt=length(time);
  t=[0;sort([time;time]);T];
  
  [rpp,cpp]=size(t);
  losses=zeros(rpp,cpp);
  
  a=parlambda(1);
  b=parlambda(2);
  switch lambda
      case 0 
      c=parlambda(3);
      ct=(1+theta)*(a*t+b/(2*pi)*cos(2*pi*c)-b/(2*pi)*cos(2*pi*(t+c)));
      case 1
      ct=(1+theta)*(a*t+b/2*t.^2);
      case 2
      c=parlambda(3);
      ct=(1+theta)*(a*t+b*(0.5*t-1/(8*pi)*sin(4*pi*(t+c))+1/(8*pi)*sin(4*pi*c)));
  end
  switch distrib
      case 'Burr'
    ct=ct*params(2)^(1/params(3))*exp(gamma(1+1/params(3))+gamma(params(1)-1/params(3))-gamma(params(1)));
      case 'exponential'
    ct=ct/params(1);
      case 'gamma'
    ct=ct*params(1)/params(2);
      case 'lognormal'
    ct=ct*exp(params(1)+params(2)^2/2);
      case 'mixofexps'
    ct=ct*(params(1)/params(2)+(1-params(1))/params(3));
      case 'Pareto'
    ct=ct*params(2)/(params(1)-1);
      case 'Weibull'
    ct=ct*exp(gamma(1+1/params(2)))/params(1)^(1/params(2));
  end
    
  
  if(lt>0)
     laux=cumsum(val);
     losses(1:2*lt+1,1)=sort([0;0;laux;laux(1:end-1)]);
     losses(2*lt+2:end)=laux(end);
  end
     
  y=u+ct-losses;
  
 ind=find(y<0,1,'first');
 if ~isempty(ind)
    y(ind:end,1)=0;
 end
