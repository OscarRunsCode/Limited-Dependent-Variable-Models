clear
clc

%Create Raq Variables
load('petdata1.mat')
train(1,:)=[];
type=table2array(train(:,1));
Age=table2array(train(:,3));
Breed1=table2array(train(:,4));
Breed2=table2array(train(:,5));
Gender=table2array(train(:,6));
Color1=table2array(train(:,7));
Color2=table2array(train(:,8));
Color3=table2array(train(:,9));
MSize=table2array(train(:,10));
FLength=table2array(train(:,11));
Vacc=table2array(train(:,12));
Worm=table2array(train(:,13));
Sterile=table2array(train(:,14));
Health=table2array(train(:,15));
Qty=table2array(train(:,16));
Fee=table2array(train(:,17));
VidQ=table2array(train(:,20));
PhotoQ=table2array(train(:,23));
SpeedR=table2array(train(:,24));

%---Create More Specific Variables---
%Health indicators
Worm2=Worm==1;       %Dewormed
Vacc2=Vacc==1;       %Vaccinated
Ster=Sterile==1;    %Sterile
Healthy=Health==1;  %Healthy
SInjr=Health==3;    %Serious Injury

%Breed Variables
%Got values from breed dictionary from dataset
Dog=Breed1<=240|Breed1==307;
Cat=Breed1>=241&Breed1<=306;

%Gender Variables
Male=Gender==1;
Female=Gender==2;
MixGen=Gender==3;  %Mixed genders like when multiple were adopted

SevAdpt=Qty>1;  %Several adoptions simultaneously
%Note that the sum of MGender=/=the sum of MixGen. This suggests problems
%in data collection. 

%Adoption Speed
SpdDy=SpeedR==0;    %Adopted 1st day
SpdWk=SpeedR==1;    %Pet was adopted between 1 and 7 days (1st week) after being listed. 
SpdMt=SpeedR==2;    %Pet was adopted between 8 and 30 days (1st month) after being listed. 
SpdQt=SpeedR==3;    %Pet was adopted between 31 and 90 days (1st quarter) after being listed. 
SpdNo=SpeedR==4;    %No adoption after 100 days of being listed. (There are no pets in this dataset that waited between 90 and 100 days).
Adopt=SpdDy+SpdWk+SpdMt+SpdQt;  %Pet was adopted


%Size at maturity
Small=MSize==1;
Medium=MSize==2;
Large=MSize==3;
XLarge=MSize==4;
NaN=MSize==0;

%Fur Size
Short=FLength==1;
MedFur=FLength==2;
LFur=FLength==3;
FurMis=FLength==0;

%---Visual---
%Summing the number adopted to later create plot
DySum=sum(SpdDy);
WkSum=sum(SpdWk);
MtSum=sum(SpdMt);
QtSum=sum(SpdQt);
QNoAdpt=sum(SpdNo);
plotx=[0 1 2 3 4];
ploty=[DySum WkSum MtSum QtSum QNoAdpt];
bar(plotx, ploty)

%Validating
TotAdopt=sum(Adopt);
QNoAdpt+TotAdopt %Should equal 14993 (ie the size of our vectors).
max(Adopt) %Should not be greater than 1 (ie was only adopted once)

save Pet2.mat

%Beginning with model
clc
clear

%Ever adopted model

%Probit
fprintf('---Probit Model---'); 
fprintf('\n'); 
data='Pet2'; dep='Adopt'; ind='Age      Dog      Male     XLarge   Medium   Large    Short    MedFur   Vacc2    Worm2    Ster     Healthy  SInjr    Qty      Fee      VidQ     PhotoQ   ';
beta=zeros(18,1);
[betaP, vc1]=probit(data,dep,ind,beta);

%Logit
fprintf('\n'); 
fprintf('---Logit Model---'); 
fprintf('\n'); 
data='Pet2'; dep='Adopt'; ind='Age      Dog      Male     XLarge   Medium   Large    Short    MedFur   Vacc2    Worm2    Ster     Healthy  SInjr    Qty      Fee      VidQ     PhotoQ   ';
betaL=zeros(18,1);
[beta, vc2]=zlogit(data,dep,ind,beta);

load Pet2
CC=[Adopt Age Dog Male XLarge Medium Large Short MedFur Vacc2 Worm2 Ster Healthy SInjr Qty Fee VidQ PhotoQ];
CM=corr(CC);

%Probit Predict
fprintf('\n'); 
fprintf('---Probit Marginal Effects---'); 
fprintf('\n'); 
ym=probit_marg(data,dep,ind,betaP);


%Percent Correct
ZC=ym(1,1)/sum(ym(:,1));
OC=ym(2,2)/sum(ym(:,2));
Zero=100*[ZC (sum(SpdNo))/14993];
ONE=100*[OC 1-((sum(SpdNo))/14993)];
fprintf('--------------------------------------\n');
fprintf('   %% Correctly Predicted \t %% in Data \n');
fprintf('0: \t\t %2.4f \t\t\t %2.4f \n', Zero);
fprintf('1: \t\t %2.4f \t\t\t %2.4f \n', ONE);

%Logit Predict
fprintf('\n'); 
fprintf('---Logit Marginal Effects---'); 
fprintf('\n'); 
ymL=logit_marg(data,dep,ind,betaP);


%Percent Correct
LZC=ym(1,1)/sum(ymL(:,1));
LOC=ym(2,2)/sum(ymL(:,2));
LZero=100*[LZC (sum(SpdNo))/14993];
LONE=100*[LOC 1-((sum(SpdNo))/14993)];
fprintf('--------------------------------------\n');
fprintf('   %% Correctly Predicted \t %% in Data \n');
fprintf('0: \t\t %2.4f \t\t\t %2.4f \n', LZero);
fprintf('1: \t\t %2.4f \t\t\t %2.4f \n', LONE);




% %1st Day Adoption
% %Probit
% data='Pet2'; dep='SpdDy'; ind='Age      Dog      Male     Small     Medium   Large    Short    MedFur   Vacc2    Worm2    Ster     Healthy  SInjr    Qty      Fee      VidQ     PhotoQ   ';
% beta=zeros(18,1);
% [beta, vc1]=probit(data,dep,ind,beta);
% 
% %Logit
% data='Pet2'; dep='SpdDy'; ind='Age      Dog      Male     Small     Medium   Large    Short    MedFur   Vacc2    Worm2    Ster     Healthy  SInjr    Qty      Fee      VidQ     PhotoQ   ';
% beta=zeros(18,1);
% [beta, vc2]=zlogit(data,dep,ind,beta);
