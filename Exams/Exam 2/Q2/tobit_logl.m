% Octave function that evaluates the Probit log-likelihood function.

function logl=tobit_logl(theta);
global x y c km;                          % Identifies global variables.
beta = theta(1:km,1);                     % extract beta (coefficients) from theta
lns = theta(km+1,1);                      % get ln(s) from theta
s = exp(lns);                             % define s
xb=x*beta;                                % define x*Beta
z=(y-xb)./s;                              % Define z
xbs = xb./s;
xbs=bound(xbs);                           % bounds call to std. normal.                       
fc=normcdf(xbs);                          % Calls normal distribution function.
fcm=1-fc;
fd=normpdf(z);
logl=sum(c.*(-1.*log(s)+log(fd))+(1-c).*log(fcm));   % Computes logL.