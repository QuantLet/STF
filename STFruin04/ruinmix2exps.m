function [ y ] = ruinmix2exps(u,theta,dparameters);
   
    % u: initial capital for risk process
	% theta: security loading in insurance collective risk model
	% dparameters: matrix, composed of 2 columns containing the parameters of loss distribution, exponential parameters (first column) and weights (second column)
    
    p1=dparameters(:,1); % exponential parameters
    p2=dparameters(:,2); % weights
    p=p2(1)/p1(1)/(p2(1)/p1(1)+(1-p2(1))/p1(2));
    psii=p1(1)*(1-p)+p1(2)*p;
    r1=(psii+theta*sum(p1)-sqrt((psii+theta*sum(p1))^2-4*prod(p1)*theta*(1+theta)))/(2*(1+theta));
    r2=(psii+theta*sum(p1)+sqrt((psii+theta*sum(p1))^2-4*prod(p1)*theta*(1+theta)))/(2*(1+theta));
    y=1/((1+theta)*(r2-r1))*((psii-r1)*exp(-r1*u)+(r2-psii)*exp(-r2*u)); % ruin probability using the Laplace transform inversion
   
end

