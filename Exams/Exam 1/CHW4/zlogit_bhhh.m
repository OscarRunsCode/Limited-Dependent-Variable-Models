% Octave function that computes the BHHH update for the Probit model.
% d is the directional vector
% vc is the estimated covariance matrix
% md is the convergence criteria

function [d,vc,md] =zlogit_bhhh(beta);
global x y km;                          % Identifies global variables.
xb=x*beta;
%xb=bound(xb);                            % Calls normal density function.
mu=0;
sigma=1;
zlog=makedist('Logistic','mu',mu,'sigma',sigma); %Create logit distribution
fl=pdf(zlog,xb);                            % Calls logit density function.
fb=cdf(zlog,xb);                            % Calls logit distribution function.
fbc=1-fb; 
r1=fl./fb;
r2=fl./fbc; 
q=((y.*r1-(1-y).*r2)*ones(1,km)).*x;
sc=sum(q)';                             % Computes gradient (score).
m=q'*q;
vc=inv(m);                              % Estimated covariance matrix.
d=vc*sc;                                % BHHH directional update.
md=norm(d,'inf');                 		 % Convergence criterion.