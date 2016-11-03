function y = simNHPP(lambda,parlambda,T,N)
% SIMNHPP Non-homogeneous Poisson process.
%       Y = SIMNHPP(lambda,parlambda,T,N) generates N trajectories of the
%       non-homogeneous Poisson process with intensity specified by LAMBDA
%       (0 - sine function, 1 - linear function, 2 - sine square function)
%       with paramters in PARLAMBDA. T is the time horizon. The function
%       usues thining method.

  a=parlambda(1);
  b=parlambda(2);
  if a<=0
      error('simNHPP: a must be a positive real number');
  end
  if a+b<= 0
    error('simNHPP: b does not fulfill special condition');
  end
  if T <= 0
     error('simNHPP: T must be a positive real number');
  end
  switch lambda
      case 0 
        c=parlambda(3);
        JM=simHPP(a+b,T,N);
      case 1
        JM=simHPP(a+b*T,T,N);
      case 2
        c=parlambda(3);
        JM=simHPP(a+b,T,N);
  end
      
  [n,m]=size(JM);
  i=1;
  y=T*ones(n,N);
  while(i<=N)
    pom=JM(find(JM(:,i)<T));
    pom=pom(2:end);
    U=rand(length(pom),1);
    
    switch lambda
        case 0 
            lambdat=(a+b*sin(2*pi*(pom+c)))/(a+b);
        case 1
            lambdat=(a+b*pom)/(a+b*T);
        case 2
            lambdat=(a+b*sin(2*pi*(pom+c)).^2)/(a+b);
    end
        
    pom=pom(find(U<lambdat));
    EN=length(pom);
    y(1,i)=0;
    y(2:EN+1,i)=pom;
    i=i+1;
  end

