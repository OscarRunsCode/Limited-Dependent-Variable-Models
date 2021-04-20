% Octave function controlling estimation of a Tobit model -
% data is a string naming a .mat file dataset
% dep is a string naming the dependent variable
% ind is a string naming the regressors.   Each regressor name
% must have a length of eight characters, the last being a space.
% Ind defines binary denoting censored data

function theta=tobit5a(data,dep,bin,ind,theta); 

% ------------------------------------------------------------------------
% This section reads data.

load(data);                   % Loads data from *.mat file.
global x y c km;           % Makes data global - accessible to functions.
y=eval(['[',dep,']']);        % Selects dependent variable.
c=eval(['[',bin,']']);        % Selects binary
nobs = size(y,1);             % Determines sample size.
x=eval(['[',ind,']']);        % Selects regressors.
x = [ ones(nobs,1) x ];       % Adds intercept.
km = size(x,2);               % Determines number of regressors.
df = nobs-km;                 % Determines degrees of freedom.

% ------------------------------------------------------------------------
% Prints title to identify output.

% fprintf('Tobit Model - The dependent variable is: %s\n',dep);
% fprintf('The data set is: %s\n',data);
% fprintf('\n'); 

% ------------------------------------------------------------------------
% Initialize some matrices for BHHH iterations.

logl0=tobit_logl(theta);               % LogL at starting values.
tol=0.0001;                         % Convergence tolerance.
maxit=1000;                          % Iteration Limit.
md=1;
i=1;

% ------------------------------------------------------------------------
% Begin BHHH iterations.

     while ((md >= tol) && (i<=maxit));       % Start iteration loop.
     [d,vc,md]=tobit_bhhh(theta);                % Compute BHHH update.

% ------------------------------------------------------------------------
% Choose largest step that increases function value.

     j=0;
     deltf = -1;
          while (deltf<0) && (j<=4);         % Start step reduction loop.
          step=(0.5).^j;                     % Reduce step by half.
          gamma=theta+(step.*d);              % Try new parameter vector.
          logl=tobit_logl(gamma);               % Evaluate Tobit logL.
          deltf=logl-logl0;                  % Determine function change.
          j=j+1;                             % Update step counter.
          end;                               % End step reduction loop.
     logl0=logl;                             % Update value of logL.
     theta=gamma;                             % Update theta

% ------------------------------------------------------------------------
% Print results for ith iteration.
% 
%      fprintf('Grad: %10.4f',md);
%      fprintf('    LogL: %10.4f',logl);
%      fprintf('    Size: %7.3f\n',step);

% ------------------------------------------------------------------------

     i=i+1;                                  % Update iteration counter.
     end;                                    % End iteration loop.

% ------------------------------------------------------------------------
% Compute Statistics.

stderr = sqrt(diag(vc));                     % Compute standard errors.
t = theta./stderr;                            % Compute t-statistics.
pvt=1-fcdf(t.^2,1,df);                        % Compute p-values.

% ------------------------------------------------------------------------
% Print Output

if (i>maxit), disp('Iteration Limit Exceeded 5a.'), end;

ind = [ 'Con     ' ind 'lnsig   '];
ind=reshape(ind,8,km+1)';

% fprintf('\n'); 
% fprintf('Regressor  Coefficient\t Std. Error \t t-stat       Prob>|t|\n');
% fprintf('--------------------------------------------------------------\n');
% for h=1:km+1;
% fprintf([ind(h,:) '%12.5f  %12.5f  %12.5f  %12.5f \n'], [theta(h) stderr(h) t(h) pvt(h)]);
% end;
sigma=exp(theta(km+1));
% fprintf(['Sigma   ' '%12.5f \n'], [sigma]);
%fprintf('Variance-Covariance Matrix \n');
%vc
clear global;
