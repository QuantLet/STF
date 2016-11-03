function y = simHPP(lambda,T,N)
% SIMHPP Homogeneous Poisson process.
%       Y = SIMHPP(lambda,T,N) generates N trajectories of the
%       homogeneous Poisson process with intensity LAMBDA. T is the time
%       horizon. 

  if lambda <= 0
    error('simHPP: Lambda must be a positive real number');
  end
  if T <= 0
    error('simHPP: T must be a positive real number');
  end
  EN = poissrnd(lambda*T,N,1);
max(EN) ;
  y=T*ones(max(EN)+1,N);
  i=1;
  while(i<=N)
      y(1,i)=0;
    if EN(i)>0
      y(2:EN(i)+1,i)=sort(T*rand(EN(i),1)); 
    end 
    i=i+1;
  end

