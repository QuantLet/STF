function y=Paretopdf(x,alpha,lambda)
%PARETOPDF Pareto probability density function (pdf).
%   Y = PARETOPDF(X,ALPHA,LAMBDA) returns the pdf of the Pareto
%   distribution with parameters ALPHA, LAMBDA, evaluated at the values in X.
%   For CONTROL=0 the error message is displayed, if the parmeters are
%   negative.
%   
%   The default values for ALPHA, LAMBDA and CONTROL are 1, 1, 0, respectively.


  if nargin<3
    lambda=1;
  end
  if nargin<2
    alpha=1;
  end
  if nargin<1
      error('stats:normpdf:TooFewInputs','Input argument X is undefined.');
  end
  if lambda<=0
      error('Non-positive sigma!');
  end
  if alpha<=0
      error('Non-positive alpha!');
  end
  y=zeros(size(x));
  pos=find(x>0); 
  y(pos)=alpha*lambda^alpha./(lambda+x(pos)).^(alpha+1);
  

