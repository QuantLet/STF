function [y] = simNHPPALP(lambda,parlambda,distrib,params,T,N)
  
  if(lambda ~= 0 && lambda ~= 1 && lambda~=2)
  	error('simNHPPALP: Lambda must be either 0,1 or 2.');
  end
  if(T <= 0 || (length(T))~=1)
  	error('simNHPPALP: T must be a positive scalar.');
  end
  if(N <= 0 || (length(N))~=1)
  	error('simNHPPALP: N must be a positive scalar.');
  end
  if(length(parlambda)~=3 && (lambda)~=1)
  	error('simNHPPALP: for lambda 0 or 2, parlambda must be a 3 x 1 vector.');
  end
  if(length(parlambda)~=2 && (lambda)==1)
  	error('simNHPPALP: for lambda 1, parlambda must be a 2 x 1 vector.');
  end
  if((strcmp(distrib,'Burr') || strcmp(distrib,'mixofexps')) && (length(params)~=3))
   	error('simNHPPALP: for Burr and mixofexps distributions, params must be a 3 x 1 vector.');
  end
  if((strcmp(distrib,'gamma') || strcmp(distrib,'lognormal')|| strcmp(distrib,'Pareto') || strcmp(distrib,'Weibull')) && (length(params))~=2)
   	error('simNHPPALP: for gamma, lognormal, Pareto and Weibull distributions, params must be a 2 x 1 vector.');
  end
  if(strcmp(distrib,'exponential') && (length(params))~=1)
   	error('simNHPPALP: for exponential distribution, params must be a scalar.');
  end
  if(strcmp(distrib, 'exponential')==0 && strcmp(distrib, 'gamma')==0 && strcmp(distrib, 'mixofexps')==0 && strcmp(distrib,'Weibull')==0 && strcmp(distrib, 'lognormal')==0 && strcmp(distrib,'Pareto')==0 && strcmp(distrib,'Burr')==0)
   	error('simNHPPALP: distribs should be: exponential, gamma, mixofexps, Weibull, lognormal, Pareto or Burr');
  end
  poisproc = simNHPP(lambda,parlambda,T,N);
  rpp      = size(poisproc,1);
  cpp      = size(poisproc,2);
  losses   = zeros(rpp,cpp);
  
  switch distrib
      case 'Burr'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux            = cumsum(Burrrnd(params(1),params(2),params(3),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:(aux-2))/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i)=zeros(rpp,1);
      end
      i = i + 1;
    end
      case'exponential'
    i = 1;
    while(i<=N)
        aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux = cumsum(exprnd(1/params(1),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i)=zeros(rpp,1);
      end
      i = i + 1;
    end
      case 'gamma'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux = cumsum(gamrnd(params(1),(1/params(2)),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i) = zeros(rpp,1);
      end
      i = i + 1;
    end
      case 'lognormal'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux = cumsum(lognrnd(params(1),params(2),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i) = zeros(rpp,1);
      end
      i = i + 1;
    end
      case 'mixofexps'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux = cumsum(mixexprnd(params(1),params(2),params(3),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i) = zeros(rpp,1);
      end
      i = i + 1;
    end
      case 'Pareto'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux = cumsum(Paretornd(params(1),params(2),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i) = zeros(rpp,1);
      end
      i = i + 1;
    end
      case 'Weibull'
    i = 1;
    while(i<=N)
    	aux = min(find(poisproc(:,i,1)==T,i,'first'));
      if(aux>2)
        laux=cumsum(wblrnd(params(1)^(-1/params(2)),params(2),aux/2-1,1));
        losses(3:aux,i) = laux(ceil((1:aux-2)/2));
        if(aux<rpp)
          losses((aux+1):rpp,i) = laux(length(laux))*ones(rpp-aux,1);
        end
      else
        losses(:,i) = zeros(rpp,1);
      end
      i = i + 1;
    end
  end
   	y        = zeros(size(poisproc));
  	y(:,:,1) = poisproc(:,:,1);
    y(:,:,2) = losses;
  
end
