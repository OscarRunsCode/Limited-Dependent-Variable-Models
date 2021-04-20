% Octave function that evaluates the Weibull log-likelihood function.

function logl=weibull_logl(theta);
global x y c km;                          % Identifies global variables.
gamma=theta(1:km,1);                       % Extracts theta from beta
beta=theta(km+1,1);                         % Extracts ln(s) from beta
xg=x*gamma;
ly=log(y);
lyb=beta.*ly;

w=-exp(xg+lyb);

xd=(1-w).*exp(xg+lyb).*beta.*y.^(beta-1);

logl=sum(c.*(log(beta)+(beta-1).*ly+xg)+w);
