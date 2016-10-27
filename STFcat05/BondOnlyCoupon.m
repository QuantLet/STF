
function [y] = BondOnlyCoupon(C,D,T,r,lambda,parlambda,distr,params,Tmax,N)


  if(lambda ~= 0 && lambda ~= 1 && lambda~=2)
  	error('BondOnlyCoupon: Lambda must be either 0,1 or 2.');
  end
  if(length(C) ~=1)
  	error('BondOnlyCoupon: coupon payments C needs to be a scalar');
  end
  if(length(r) ~=1)
  	error('BondOnlyCoupon: discount rate needs to be a scalar');
  end
  if(length(D)==1)
  	error('BondOnlyCoupon: threshold level D needs to be a vector ');
  end
  if(length(T)==1)
  	error('BondOnlyCoupon: time to expiry T needs to be a vector ');
  end
  x  = simNHPPALP(lambda,parlambda,distr,params,Tmax,N);
  Tl = length(T);
  Dl = length(D);
  y  = zeros(Tl*Dl,3);
  i   = 1; %loop (times to maturity)
  j   = 1; %loop (treshold levels)
  k   = 1; %loop (trajectories)
  wyn = 0;
  while(i<=Tl)
    while(j<=Dl)
      while(k<=N)
        traj  = [x(:,k,1),x(:,k,2)];
        if (traj(length(traj(find(traj(:,1)<=T(i)),1)),2)<=D(j))
          wyn = wyn + (1-exp(-r*T(i)))/r;
        else
          wyn = wyn + (1-exp(-r*traj(length(traj(find(traj(:,2)<=D(j)))),1)))/r;
        end
        k = k + 1;
      end
      y((i-1)*Dl+j,1) = T(i);
      y((i-1)*Dl+j,2) = D(j);
      y((i-1)*Dl+j,3) = C*wyn/N;
      wyn = 0;
      k   = 1;
      j   = j + 1;
    end
    j = 1;
    i = i + 1;
  end

end
