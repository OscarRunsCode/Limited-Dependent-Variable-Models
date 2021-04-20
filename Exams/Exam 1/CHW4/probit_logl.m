% Octave function that evaluates the Probit log-likelihood function.

function logl=probit_logl(beta);
global x y km;                          % Identifies global variables.
xb=x*beta;
%xb=bound(xb);                            % Calls normal distribution function.
fb=normcdf(xb);                            % Calls normal distribution function.
fbc=1-fb; 
logl=sum(y.*log(fb)+(1-y).*log(fbc));   % Computes logL.