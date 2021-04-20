% Octave function that computes the BHHH update for the Weibull model.
% d is the directional vector
% vc is the estimated covariance matrix
% md is the convergence criteria

function [d,vc,md] =weibull_bhhh(theta);
global x y c km;                          % Identifies global variables.
gamma=theta(1:km,1);                       % Extracts theta from beta
beta=theta(km+1,1);                         % Extracts ln(s) from beta
xg=x*gamma;
ly=log(y);
lyb=beta.*ly;

qgamma=((c-exp(xg+lyb))*ones(1,km)).*x;
qbeta=c./beta+(c-exp(xg+lyb)).*ly;
q=[qgamma qbeta];

sc=sum(q)';
m=q'*q;
vc=inv(m);                              % Estimated covariance matrix.
d=vc*sc;                                % BHHH directional update.
md=norm(d,'inf');                 		 % Convergence criterion.