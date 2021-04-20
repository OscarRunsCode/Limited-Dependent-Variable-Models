% Octave function that estimates marginal effects in a Probit model -
% data is a string naming a .mat file dataset
% dep is a string naming the dependent variable
% ind is a string naming the regressors.   Each regressor name
% must have a length of eight characters, the last being a space.

% This program provides point estimate and standard errors for f(Xd)d,
% where f( ) is the standard normal density function, X is the vector of
% sample means of the regressors, and d is the vector of Probit 
% coefficient estimates.  The reported "marginal effect" of the intercept
% is obviously nonsence.  The reported "marginal effect" of any binary
% regressors should be a useful approximation of the "incremental effect"
% evaluated at the means of the regressors.

function probit_marg(data,dep,ind,beta); 
clc

% ------------------------------------------------------------------------
% This section reads data.

load(data);                   % Loads data from *.mat file.
global x y km;                % Makes data global - accessible to functions.
y=eval(['[',dep,']']);        % Selects dependent variable.
nobs = size(y,1);             % Determines sample size.
y = (y>0);                    % Insures that dependent variable is binary.
x=eval(['[',ind,']']);        % Selects regressors.

x = [ ones(nobs,1) x ];       % Adds intercept.
km = size(x,2);               % Determines number of regressors.
df = nobs-km;                 % Determines degrees of freedom.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('Probit Marginal Effects - The dependent variable is: %s\n',dep);
fprintf('The data set is: %s\n',data);
fprintf('\n'); 

% ------------------------------------------------------------------------
% Compute marginal effects at sample means.

xbar=sum(x)./nobs;					  % Compute sample means of regressors.
xb=x*beta;
xbb=xbar*beta;
%xb=bound(xb);
%xbb=bound(abb);
fl=normpdf(xbb);                       % Calls normal density function.
fb=normcdf(xb);                       % Calls normal distribution function.
marg=fl.*beta;

% ------------------------------------------------------------------------
% Compute covariance matrix of marginal effects.

[d,vc,md]=probit_bhhh(beta);           % Compute VC matrix of coefficients.

delt=eye(km)-xbb.*(beta*xbar);
delt=fl.*delt;
vc=delt*vc*(delt');

% ------------------------------------------------------------------------
% Compare predicted and actual values.

yp=fb>0.5;                     % Predicted y=1.
ya=[ 1-y y ];
yb=[ 1-yp yp ];
ym=ya'*yb;                     % Prediction matrix.

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = marg./stderr;                            % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df);                        % Compute p-values.

% ------------------------------------------------------------------------
% Print Output

ind = [ 'Con     ' ind];
ind=reshape(ind,8,km)';

fprintf('\n'); 
fprintf('Regressor   Marginal\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:km;
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [marg(h) stderr(h) t(h) pvt(h)]);
end;

fprintf('\n'); 
fprintf('\n'); 
fprintf('\t \t Predicted 0\t Predicted 1\n');
fprintf(['Actual 0\t %7.0f \t  %7.0f \n'], ym(1,:));
fprintf(['Actual 1\t %7.0f \t  %7.0f \n'], ym(2,:));


clear global;
endfunction;