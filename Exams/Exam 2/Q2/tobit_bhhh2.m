% Octave function that computes the BHHH update for the Tobit model.
% d is the directional vector
% vc is the estimated covariance matrix
% md is the convergence criteria

function [d,vc,md] =tobit_bhhh2(theta);
global x y c km;                          % Identifies global variables.
beta = theta(1:km,1);
s  = theta(km+1,1);
%s = exp(lns);
% Define Tobit Score equation components
xb=x*beta;                                  % Calls normal density function.
xbs=xb./s;
xbs=bound(xbs);                               % bounds call to std. normal.
fl=normpdf(xbs);                            % Calls normal density function.
fb=normcdf(xbs);                            % Calls normal distribution function.
fbc=1-fb;
z =(y-xb)./s;
e1=z./s;
e2=(fl./fbc)./s;                        %Inv Mills ratio * s^-1
qb=((c.*e1-(1-c).*e2)*ones(1,km)).*x;   %First set of score equations (for beta)
e3=-1./s+(z.*z)./s;                     %First term in score eqn for s (without c)
e4=xbs.*(fl./fbc)./s;                   %Second term in score eqn for s (without 1-c)
qs=(c.*e3+(1-c).*e4);                   % Score eqn for sigma
q=[qb qs];                              % Stack k+1 score equations;
sc=sum(q)';                             % Computes gradient (score).
m=q'*q;
vc=inv(m);                              % Estimated covariance matrix.
d=vc*sc;                                % BHHH directional update.
md=norm(d,'inf');                 		 % Convergence criterion.