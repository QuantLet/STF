function y=mixexppdf(x,alpha,beta1,beta2,control)
%MIXEXPPDF Mixed exponential probability density function (pdf).
%   Y = MIEXEXPPDF(X,ALPHA,BETA1,BETA2) returns the pdf of the mixed
%   exponential distribution with mixing probability A and distributions 
%   parameters BETA1, BETA2, evaluated at the values in X.
%   For CONTROL=0 the error message is displayed, if the parmeters are
%   negative or a>1.
%   
%   The default values for A, BETA1, BETA2 and CONTROL are 0.5, 1, 2, 0
%   respectively.


  if nargin<5
      control=0;
  end
  if nargin<2
      a=.5;
  end
  if nargin<3
    b1=1;
  end
  if nargin<4
    b2=2;
  end
  if(control==0)
      if beta1<=0
        error('Non-positive beta1!');
      end
      if beta2<=0
        error('Non-positive beta2!');
      end  
      if alpha<=0
          error('Alpha lesser or equal 0!');
      end
      if alpha>=1
          error('Alpha greater or equal 1!');
      end
  end

  
  y=zeros(size(x));
  pos=find(x>0);
  y(pos)=alpha*beta1*exp(-beta1*x(pos))+(1-alpha)*beta2*exp(-beta2*x(pos));
