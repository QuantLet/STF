function [x0,y0] = empcdf(x,infsupport)
%EMPCDF Empirical cumulative distribution function (cdf).
%   EMPCDF(X) plots the empirical cdf of the elements in vector X assuming 
%   that the support of the distribution is (-INF,MAX(X)].
%   EMPCDF(X,INFSUPPORT) allows the user to specify the infimum of the 
%   support.
%   [X0,Y0] = EMPCDF(X) does not draw a graph, but returns vectors X0 and 
%   Y0 such that PLOT(X0,Y0) is the empirical cdf.
%
%   See also EMPDENS.

%   Written by Rafal Weron (1999.09.15, rev. 2004.04.03)
%   Copyright (c) 1999-2006 by Rafal Weron

if nargin == 0
    error('Requires one input argument.')
end

if nargin < 2
    infsupport = -inf;
end

x = x(:);
n = length(x);
x0 = sort([infsupport; x]);
y0 = (0:n)'/n;

if nargout == 0
    plot(x0,y0);
end