% Octave function for Heckman-Lee estimates of a Selection model -
% data is a string naming a .mat file dataset
% dep is a string naming the dependent variable
% bin is a string naming the indicator variable
% ind1 is a string naming the regressors of the selection equation.
% ind2 is a string naming the regressors of the regression equation.
% Each regressor name  must have a length of eight characters, the last being a space.

function [t0,t1,t2]=select_hl(data,bin,ind1,dep,ind2); 
clc

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('Heckman-Lee Estimates of Selection Model\n');
fprintf('The data set is: %s\n',data);
fprintf('\n'); 


% ------------------------------------------------------------------------
% Load data for second stage

load(data);                   % Loads data from *.mat file.
c=eval(['[',bin,']']);        % Selects indicator variable.
n1=sum(c);          	      % Determines subsample size.
y=eval(['[',dep,']']);        % Selects dependent variable.
nobs = size(y,1);             % Determines sample size.

w=eval(['[',ind1,']']);        % Selects regressors for selection equation.
w = [ ones(nobs,1) w ];       % Adds intercept.
kw = size(w,2);                % Determines number of regressors.

x=eval(['[',ind2,']']);        % Selects regressors for regression equation.
x = [ ones(nobs,1) x ];       % Adds intercept.
kx = size(x,2);                % Determines number of regressors.
df = n1-kx;                 % Determines degrees of freedom.

% ------------------------------------------------------------------------
% First-Stage Probit Estimates

a=zeros(kw,1);
a=probit(data,bin,ind1,a);

% ------------------------------------------------------------------------
% Construct auxilary regressor for second stage.

wa=w*a;
%wa=bound(wa);
fl=normpdf(wa); 
fb=normcdf(wa);
r1=fl./fb;

% ------------------------------------------------------------------------
% Determine observed subsample.

y1=y(c>0);
x1=x(c>0,:);
r1=r1(c>0);
wa=wa(c>0);

% ------------------------------------------------------------------------
% Probit-OLS-Rho=0.

b=inv(x1'*x1)*(x1'*y1);
s=sqrt((y1-x1*b)'*(y1-x1*b)/df);
vc=(s.^2).*inv(x1'*x1);
logl=sum(log(normpdf((y1-x1*b)./s)./s));

ls=log(s);
t0=[a;b;ls];
t1=[a;b;ls;0];

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = b./stderr;                           % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df);                        % Compute p-values.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('\n'); 
fprintf('OLS Coefficients - The dependent variable is: %s\n',dep);
fprintf('\n'); 

% ------------------------------------------------------------------------
% Print Output

ind = [ 'Con     ' ind2 ];
ind=reshape(ind,8,kx)';

fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:kx;
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [b(h) stderr(h) t(h) pvt(h)]);
end;
fprintf('ln(Sig) %12.5f\n',ls);
ro=0;
fprintf('FisherZ %12.5f\n',ro);
fprintf('LogL    %12.5f\n',logl);

% ------------------------------------------------------------------------
% Heckman-Lee

x1=[x1 r1];
b=inv(x1'*x1)*(x1'*y1);
s=sqrt((y1-x1*b)'*(y1-x1*b)/(df-1));
vc=(s.^2).*inv(x1'*x1);
sro=b(kx+1);
s2=(s.^2)+((sro.^2).*r1'*(wa+r1)/(df-1));
s=sqrt(s2);
ro=sro/s;

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = b./stderr;                           % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df-1);                        % Compute p-values.

% ------------------------------------------------------------------------
% Prints title to identify output.

fprintf('\n'); 
fprintf('Heckman-Lee Coefficients - The dependent variable is: %s\n',dep);
fprintf('Uncorrected Standard Errors\n');
fprintf('\n'); 

% ------------------------------------------------------------------------
% Print Output

ind = [ 'Con     ' ind2 'Mills   ' ];
ind=reshape(ind,8,(kx+1))';

fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
fprintf('--------------------------------------------------------------\n');
for h=1:(kx+1);
fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [b(h) stderr(h) t(h) pvt(h)]);
end;

% ------------------------------------------------------------------------
% Bound the HL estimate of sigma.

if s<=0;
	s=1
   fprintf('Estimate of Sigma is non-positive.\n');
end;

ls=log(s);
fprintf('Ln(Sig) %12.5f\n',ls);

% ------------------------------------------------------------------------
% Bound the HL estimate of rho.

if abs(ro)>=1;
	ro=sign(ro).*(0.9);
   fprintf('Estimate of Rho is outside the unit interval.\n');
end;

% ------------------------------------------------------------------------
% Convert Rho to Fisher-Z.

ro=log(1+ro)-log(1-ro);
fprintf('FisherZ %12.5f\n',ro);

t2=[a;b(1:kx);ls;ro];
endfunction;
