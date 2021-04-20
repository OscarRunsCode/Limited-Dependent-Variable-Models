% Octave function that evaluates the Probit log-likelihood function.

function logl=zlogit_logl(beta);
global x y km;                          % Identifies global variables.
mu=0;
sigma=1;
zlog=makedist('Logistic','mu',mu,'sigma',sigma); %Create logit distribution
xb=x*beta;
%xb=bound(xb);                            % Calls normal distribution function.
fb=cdf(zlog,xb);                            % Calls normal distribution function.
fbc=1-fb; 
logl=sum(y.*log(fb)+(1-y).*log(fbc));   % Computes logL.