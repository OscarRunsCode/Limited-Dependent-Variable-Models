% Octave function for OLS estimation -
% data is a string naming a .mat file dataset
% dep is a string naming the dependent variable
% ind is a string naming the regressors.   Each regressor name
% must have a length of eight characters, the last being a space.

function beta=ols(data,dep,ind);
clc

% ------------------------------------------------------------------------
% This section reads data.

load(data);                   % Loads data from *.mat file.
y=eval(['[',dep,']']);        % Selects dependent variable.
nobs = size(y,1);             % Determines sample size.
x=eval(['[',ind,']']);        % Selects regressors.
x = [ ones(nobs,1) x ];       % Adds intercept.
km = size(x,2);               % Determines number of regressors.
df = nobs-km;                 % Determines degrees of freedom.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('OLS Estimates - The dependent variable is: %s\n',dep);
fprintf('The data set is: %s\n',data);
fprintf('\n'); 

% ------------------------------------------------------------------------
% Obtain point estimates.

vc=inv(x'*x);
beta=vc*(x'*y);
s2=(y-x*beta)'*(y-x*beta)/df;
vc=s2.*vc;

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = beta./stderr;                            % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df);                        % Compute p-values.

% ------------------------------------------------------------------------
% Print Output

ind = [ 'Con     ' ind];
ind=reshape(ind,8,km)';

fprintf('\n'); 
fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:km;
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [beta(h) stderr(h) t(h) pvt(h)]);
end;
