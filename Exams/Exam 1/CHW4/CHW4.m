%-----Starting from scratch-----
clc
clear

%Diary
diary CHW4_Output_Oscar_Martinez.txt

%Introduction
fprintf('--------------------------------------------------------------\n');
fprintf('Oscar Martinez \t Computer Homework 4 \t Metrics III\n');
fprintf('--------------------------------------------------------------\n');

fprintf('--------------------------------------------------------------\n');
fprintf('\t \t  Probit Model \n');
fprintf('--------------------------------------------------------------\n');

%Importing the data:
load mls.txt;

%Creating Variables from columns of imported data
t = mat2cell(mls, size(mls,1), ones(1,size(mls,2)));
[sp age lot sqft beds gar mfi pmin paved fin vac trav days ap] = deal(t{:});

%Removing extraneous variables
clear t mls;

%saving as .mat
save mls.mat

%Creating new, needed variables
sold=sp>0; sqft2=sqft.^2; lsp = log(sp+(1-sold));

%Saving updated data
save mls2.mat

%Clearing workspace for models
clear

%Note: commented out the clc from the function
%-----Probit Model-----
%Creating needed variables:
data='mls2'; dep='sold'; ind='age     lot     sqft    beds    gar     mfi     pmin    paved   fin     vac     trav    ap      ';
beta=zeros(13,1);
beta=probit(data,dep,ind,beta);

%Modifying the output for ease of viewing
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\t \t  Logit Model \n');
fprintf('--------------------------------------------------------------\n');

%-----Logit Model-----
%Resetting the beta
beta=zeros(13,1);

%Running the function
beta=zlogit(data,dep,ind,beta);

%Aesthetics
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf('--------------------------------------------------------------\n');

%closing output
diary off

%-----Changes made to function files-----
%The changes made, other than renaming anything the prefix "prob" to "(z)log" and commenting out the clc, were:
% -> In zlogit_bhhh.m:
    % Created variables mu and  sigma with values 0 and 1, respectively.
    % Used "zlog=makedist('Logistic','mu',mu,'sigma',sigma);" to define a Logistic distribution. 
    % Replaced "fl=normpdf(xb);" with "fl=pdf(zlog,xb);" and "fb=normcdf(xb);" with "fb=cdf(zlog,xb);".
% -> In zlogit_logl.m:
    % Created variables mu and  sigma with values 0 and 1, respectively.
    % Used "zlog=makedist('Logistic','mu',mu,'sigma',sigma);" to define a Logistic distribution. 
    % Replaced "fb=normcdf(xb);" with "fb=cdf(zlog,xb);".