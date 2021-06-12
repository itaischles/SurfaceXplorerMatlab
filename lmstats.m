function [var,std,C] = lmstats(FUN,x)
% [var,std,C] = lmstats(FUN,x)
% FUN - the anonymous function which outputs the residual at the variables
%       estimates
% x - point at which we wish to calculate the standard deviation
% 
% MSE - Mean Square Error
% std - standard deviation of x
% C - the covariance matrix
% ----------------------------
% Matthew Krzyaniak
% Northwestern University
% 2016
% ----------------------------
if nargin == 2
  weights = 1;
end

% First calculate the Mean Square Error
fit = feval(FUN,x);
var = (fit(:)'*fit(:))/(length(fit(:))-1);

% Calculate the jacobian matrix using a central difference approximation
lx = length(x);
J  = zeros(length(fit(:)),lx);
for k = 1:lx
  if x(k)==0
    dx = nthroot(eps,3);
  else
    dx = x(k)*nthroot(eps,3);
  end
  xdp = x;
  xdm = x;
  xdp(k) = xdp(k)+dx;
  xdm(k) = xdm(k)-dx;
  rd = feval(FUN,xdp);
  rc = feval(FUN,xdm);
  J(:,k)=(rd-rc)/(2*dx);
end
Hes = (J'*diag(weights(:))*J);
% Calculate the covariance matrix
C = var*inv(Hes);
% calculate the standard devation of the parameters
std = sqrt(diag(C));