% Octave function for Heckman-Lee estimates of a Tobit model -
% data is a string naming a .mat file dataset
% dep is a string naming the dependent variable
% bin is a string naming the indicator variable
% ind is a string naming the regressors.
% Each regressor name  must have a length of eight characters, the last being a space.

function [t1,t2]=tobit_hl(data,dep,bin,ind); 
clc

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('Heckman-Lee Estimates of Tobit Model\n');
fprintf('The data set is: %s\n',data);
fprintf('\n'); 


% ------------------------------------------------------------------------
% Load data for second stage

load(data);                   % Loads data from *.mat file.
c=eval(['[',bin,']']);        % Selects indicator variable.
n1=sum(c);          	      % Determines subsample size.
y=eval(['[',dep,']']);        % Selects dependent variable.
nobs = size(y,1);             % Determines sample size.
x=eval(['[',ind,']']);        % Selects regressors.

x = [ ones(nobs,1) x ];       % Adds intercept.
kx = size(x,2);                % Determines number of regressors.
df = n1-kx;                 % Determines degrees of freedom.

% ------------------------------------------------------------------------
% First-Stage Probit Estimates
a=zeros(kx,1);
a=probit(data,bin,ind,a);

% ------------------------------------------------------------------------
% Construct auxilary regressor for second stage.

xa=x*a;
%xa=bound(xa);
fl=normpdf(xa); 
fb=normcdf(xa);
r1=fl./fb;

% ------------------------------------------------------------------------
% Determine observed subsample.

y1=y(c>0);
x1=x(c>0,:);
r1=r1(c>0);
xa=xa(c>0);

% ------------------------------------------------------------------------
% OLS.

b=inv(x1'*x1)*(x1'*y1);
s=sqrt((y1-x1*b)'*(y1-x1*b)/df);
vc=(s.^2).*inv(x1'*x1);
logl=sum(log(normpdf((y1-x1*b)./s)./s));
ls=log(s);
t1=[b;ls];

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = b./stderr;                           % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df);                        % Compute p-values.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('\n'); 
fprintf('\n'); 
fprintf('OLS Coefficients - The dependent variable is: %s\n',dep);

% ------------------------------------------------------------------------
% Print Output

ind = [ 'Con     ' ind ];
ind=reshape(ind,8,kx)';

fprintf('\n'); 
fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:kx;
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [b(h) stderr(h) t(h) pvt(h)]);
end;
fprintf('Ln(Sig) %12.5f\n',ls);
fprintf('LogL    %12.5f\n',logl);

% ------------------------------------------------------------------------
% Heckman-Lee

x1=[x1 r1];
b=inv(x1'*x1)*(x1'*y1);
s=sqrt((y1-x1*b)'*(y1-x1*b)/(df-1));
vc=(s.^2).*inv(x1'*x1);
sro=b(kx+1);
s2=(s.^2)+((sro.^2).*r1'*(xa+r1)/(df-1));
s=sqrt(s2);

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = b./stderr;                           % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df-1);                        % Compute p-values.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('\n'); 
fprintf('\n'); 
fprintf('Heckman-Lee Coefficients - The dependent variable is: %s\n',dep);
fprintf('Uncorrected Standard Errors\n');

% ------------------------------------------------------------------------
% Print Output

ind = [ ind ;'Mills   ' ];

fprintf('\n'); 
fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:(kx+1);
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [b(h) stderr(h) t(h) pvt(h)]);
end;

if s<=0;
   s=1;
   fprintf('Estimate of Sigma is non-positive.\n');
end;

ls=log(s);
fprintf('Ln(Sig) %12.5f\n',ls);

t2=[b(1:kx);ls];
endfunction
