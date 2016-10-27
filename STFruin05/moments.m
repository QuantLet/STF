% returns the k-th moment (up to fourth) of the mixture of 2 exponentials distribution claims

function [ mk ] = moments( k,dparameters );

    % k: order of moment to calculate
	% dparameters: list, composed of 2 vectors containing the parameters of loss distribution, weights (first vector) and exponential parameters (second vector)
    p1=dparameters(:,1); % weights
    p2=dparameters(:,2); % exponential parameters
   if (k==1)
    mk=sum(p1./p2);
   elseif (k==2)
    mk=2*sum(p1./(p2.^2));
   elseif (k==3)
	mk=6*sum(p1./(p2.^3));
   elseif (k==4)
    mk=24*sum(p1./(p2.^4));
   else
       disp('k chosen too large')
end

