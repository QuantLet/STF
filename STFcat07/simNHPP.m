function [y] = simNHPP(lambda,parlambda,T,N)
  	
% SIMNHPP Non-homogeneous Poisson process.
%       Y = SIMNHPP(lambda,parlambda,T,N) generates N trajectories of the
%       non-homogeneous Poisson process with intensity specified by LAMBDA
%       (0 - sine function, 1 - linear function, 2 - sine square function)
%       with paramters in PARLAMBDA. T is the time horizon. The function
%       usues thining method.

  a = parlambda(1);
  b = parlambda(2);
  if (a<=0)
      error('simNHPP: a must be a positive real number');
  end
  if (a+b<= 0)
      error('simNHPP: b does not fulfill special condition');
  end
  if (T <= 0)
      error('simNHPP: T must be a positive real number');
  end
  if (lambda == 0)
        c  = parlambda(3);
        JM = simHPP(a+b,T,N);
  elseif( lambda==1)
        JM = simHPP(a+b*T,T,N);
  elseif(lambda==2)
        c  = parlambda(3);
        JM = simHPP(a+b,T,N);
  end
  
  rjm      = size(JM,1);
  y        = ones(rjm,N,2);
  y(:,:,1) = T*y(:,:,1);
  y(:,:,2) = 0*y(:,:,2);
  maxEN    = 0;
  i        = 1;
  while(i<=N)
  	JMind = find(JM(:,i,1)<T);
    pom   = JM(JMind,i,1);
    poml  = length(pom);
    pom   = pom(2*(1:(poml-1)/2));
    U     = unifrnd(0,1,length(pom),1);
    
    if(lambda == 0)
            lambdat = (a+b*sin(2*pi*(pom+c)))/(a+b);
    elseif(lambda == 1)
            lambdat = (a+b*pom)/(a+b*T);
    elseif(lambda == 2)
            lambdat = (a+b*sin(2*pi*(pom+c)).^2)/(a+b);
    end
    pomind            = find(U<lambdat);
    pom               = pom(pomind);
    EN                = length(pom);
    maxEN             = max([maxEN,EN]);
    y(1:(2*EN+1),i,1) = [0;pom(ceil((1:(2*EN))/2))];
    y(2:(2*EN),i,2)   = floor((1:(2*EN-1))/2)';
    y((2*EN+1):rjm,i,2) = EN*ones((rjm-2*EN),1);
    i               = i+1;
  end
  y = y(1:(2*maxEN+2),:,:);
  
end