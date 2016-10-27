% returns the moment generating function or its k-th derivative (up to third) for mixture of 2 exponentials distribution claims

function [ y ] = mgfs( x,k,dparameters )

	% x: scalar, n x 1 vector or m x n matrix,  argument of the moment generating function
	% k: scalar, integer, 0 =< k <= 3, order of the derivative
	% dparameters: vector, composed of 2 columns containing the parameters of the loss distribution, weights (first column) and exponential parameters (second column)
	
    p1=dparameters(:,1); % weights
    p2=dparameters(:,2); % exponential parameters
   if (k==0)
    y=sum((p1.*p2)./(p2-x'));
   elseif(k==1)
    y=sum((p1.*p2)./((p2-x').^2));
   elseif(k==2)
	y=2*sum((p1.*p2)./((p2-x').^3));
   elseif(k==3)
	y=6*sum((p1.*p2)./((p2-x').^4));
end

