function [t,y]=simNHPPRP(u,theta,lambda,parlambda,distrib,params,T,N)
% SIMNHPPRP Simulated trajectories of the risk process.
%       [t,y]=SIMNHPPRP(u,theta,lambda,parlambda,distrib,params,T,N) 
%       returns N simulated trajectories of the risk process with claim 
%       sizes from the distribution specified in DISTRIB with parameters 
%       in PARAMS. The claim arrival process corresponds to the
%       non-homogeneous Poisson process with inensity specified by LAMBDA 
%       (0 - sine function, 1 - linear function, 3 - sine square function) 
%       with parameters in PARLAMBDA. THETA is the relative safety loading,
%       U the initial capital and T the time horizon.

  poisproc=simNHPP(lambda,parlambda,T,N);
  [rpp,cpp]=size(poisproc);
  losses=zeros(2*rpp-1,cpp);
  
  a=parlambda(1);
  b=parlambda(2);
  
switch lambda
      case 0 
      c=parlambda(3);
      ct=(1+theta)*(a*poisproc+b/(2*pi)*cos(2*pi*c)-b/(2*pi)*cos(2*pi*(poisproc+c)));
      case 1
      ct=(1+theta)*(a*poisproc+b/2*poisproc.^2);
      case 2
      c=parlambda(3);
      ct=(1+theta)*(a*poisproc+b*(0.5*poisproc-1/(8*pi)*sin(4*pi*(poisproc+c))+1/(8*pi)*sin(4*pi*c)));
end
  
  switch distrib
      case 'Burr'
    ct=ct*params(2)^(1/params(3))*exp(gamma(1+1/params(3))+gamma(params(1)-1/params(3))-gamma(params(1)));
    i=1;
        while(i<=N)
          aux=find(poisproc(:,i)==T,1,'first');
            if(aux>0)
                laux=cumsum(Burrrnd(params(1),params(2),params(3),aux-1,1));
                 losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
                 losses(2*aux:end,i)=laux(end);
            else
                losses(1:rpp,i)=zeros(rpp,1);
            end             
      i=i+1;
    end
      case 'exponential' 
    ct=ct/params(1);
     i=1;
    while(i<=N)
      aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(exprnd(1/params(1),aux-1,1));
            losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
            losses(2*aux:end,i)=laux(end);
            end            
      i=i+1;
    end
      case 'gamma'
    ct=ct*params(1)/params(2);
     i=1;
    while(i<=N)
      aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(gamrnd(params(1),1/params(2),aux-1,1));
          losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
          losses(2*aux:end,i)=laux(end);
            end            
      i=i+1;
    end
      case 'lognormal'
    ct=ct*exp(params(1)+params(2)^2/2);
   i=1;
    while(i<=N)
      aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(lognrnd(params(1),params(2),aux-1,1));
         losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
            losses(2*aux:end,i)=laux(end);
            end            
      i=i+1;
    end
      case 'mixofexps'
    ct=ct*(params(1)/params(2)+(1-params(1))/params(3));
    i=1;
    while(i<=N)
     aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(mixexprnd(params(1),params(2),params(3),aux-1,1));
         losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
           losses(2*aux:end,i)=laux(end);
            end               
      i=i+1;
    end
      case 'Pareto'
   ct= ct*params(2)/(params(1)-1);
     i=1;
    while(i<=N)
      aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(Paretornd(params(1),params(2),aux-1,1));
          losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
           losses(2*aux:end,i)=laux(end);
            end             
      i=i+1;
    end
      case 'Weibull'
    ct=ct*exp(gamma(1+1/params(2)))/params(1)^(1/params(2));
    i=1;
    while(i<=N)
     aux=find(poisproc(:,i)==T,1,'first');
      if(aux>0)
        laux=cumsum(wblrnd(params(1)^(-1/params(2)),params(2),aux-1,1));
         losses(1:2*aux-1,i)=sort([0;0;laux;laux(1:end-1)]);
          losses(2*aux:end,i)=laux(end);
            end             
      i=i+1;
    end
  end

  y=u+sort([ct(2:end,:);ct])-losses;
  
  for j=1:N
      ind=find(y(:,j)<0,1,'first');
      if ~isempty(ind)
          y(ind:end,j)=0;
      end
  end
  t=sort([poisproc(2:end,:);poisproc(:,:)]);

