 % Computer Take-Home 2: Question 2
 % Metrics III
 %Oscar Martinez
 
 %Diary
 diary Q1_Output_Oscar_Martinez.txt
 
 %Introduction
 fprintf('--------------------------------------------------------------\n');
 fprintf('Oscar Martinez \t Take-Home 2: Question 2 \t Metrics III\n');
 fprintf('--------------------------------------------------------------\n');
 
 clear
 clc
 load martins;
 % shorten length of pexpchd2;
 pexpch2=pexpchd2;
 save martins1;
 clear;
 
 % Now identify common threshold and transform lwage
 load martins1;
 nemploy=(employ==0);
 nobs=size(employ,1);
 
 % create vector of maximums ;
 mins=ones(nobs,1).*max(lwage);
 lwage2=lwage.*employ+mins.*nemploy;
 % generate lambda;
 lambda=min(lwage2)
 % create transformed dependent variable as lwage2;
 % lwage2 equals lwage minus the minimum of non-zero ;
 % values when employ=1 and 0 otherwise;
 lwage2=(lwage2-lambda).*employ;
 wage2=exp(lwage2);
 clear lambda nobs ans;
 save martins2;
 clear;
 
 %define text strings for tobit models;
 data='martins2';
 dep2='lwage2';
 ind='edu     pexp    pexp2   pexpchd pexpch2 ';
 bin='employ';
 
 % Use HL to generate startingvalues;
 [s1b,s2b]=tobit_hl(data,dep2,bin,ind);
 theta=tobit2a(data,dep2,bin,ind,s2b);
 theta(7)=exp(theta(7));
 fprintf('Common Threshold Model (based on Transformed lwage)');
 [theta, vc,stderr]=tobit3(data,dep2,bin,ind,theta);
 
 %95% Confidence Interval
 BCo=theta(2);
 BErr=stderr(2);
 LBCI=BCo-1.96*BErr;   %Lower CI
 HBCI=BCo+1.96*BErr;   %Upper CI
 fprintf('\n');
 fprintf('\n-----------Part A--------------- \n');
 fprintf('---Confidence Interval---\n');
 fprintf('(%2.4f,%2.4f)',LBCI, HBCI);
 
 %Part B
 fprintf('\n');
 fprintf('\n-----------Part B--------------- \n');
 fprintf('---Variance-Covariance Matrix---\n');
 vc
 fprintf('    Alpha      Beta     Delta     Gamma     Theta     Lambda     Sigma \n');
 
 % get stats needed for marginal effect
 load martins
 lwage(lwage==0)=NaN;
 M=[child, pexp, edu, lwage];
 fprintf('\n---Median Values---');
 MED=nanmedian(M);
 MWAGE=exp(MED(1,4)); MPEXP=MED(1,2); MCHD=MED(1,1); MEDU=MED(1,3);
 MEDout=[MED MWAGE]
 fprintf('    child     pexp        edu     lwage     wage\n');
 
 
 fprintf('\n---Marginal Effect---');
 %Marginal Effect
 M1=theta(3)+2*theta(4)*MPEXP+theta(5)*MCHD+2*theta(6)*MPEXP*MCHD; %Marginal Effect w/o lwage
 MargEff=M1*MWAGE
 
 %Variances
 fprintf('\n---Std Err---');
 VB=vc(2,2); VD=vc(3,3); VG=vc(4,4); VT=vc(5,5); VL=vc(6,6); CDG=vc(3,4); CDT=vc(3,5);
 CDL=vc(3,6); CGT=vc(4,5); CGL=vc(4,6); CTL=vc(5,6);
 
 %Could've combined them all into one equation but this is easier to error
 %correct
 VARR1=(VD+4*VG*MPEXP^2+MCHD^2*VT+4*MPEXP^2*MCHD^2*VL);
 VARR2=2*2*MPEXP*CDG+2*MCHD*CDT+2*2*MPEXP*MCHD*CDL;
 VARR3=2*2*MPEXP*MCHD*CGT+2*2*MPEXP*2*MPEXP*MCHD*CGL;
 VARR4=2*MCHD*2*MPEXP*MCHD*CTL;
 VSUM=VARR1+VARR2+VARR3+VARR4;
 VMARG=MWAGE^2*VSUM;
 MargSD=sqrt(VMARG)
 
 %95% Confidence Interval
 LBCI2=MargEff-1.96*MargSD;   %Lower CI
 HBCI2=MargEff+1.96*MargSD;   %Upper CI
 fprintf('---Confidence Interval---\n');
 fprintf('(%2.4f,%2.4f)',LBCI2, HBCI2);
 
 
 fprintf('\n');
 fprintf('\n-----------Part C--------------- \n');
 
 %Backup Relevant Files
 theta1=theta;
 vc1=vc;
 
 %define text strings for tobit models;
 data='martins2';
 dep2='wage2';
 ind='edu     pexp    pexp2   pexpchd pexpch2 ';
 bin='employ';
 
 % Use HL to generate startingvalues;
 [s1b,s2b]=tobit_hl2(data,dep2,bin,ind);
 theta=tobit5a(data,dep2,bin,ind,s2b);
 theta(7)=exp(theta(7));
 %fprintf('Common Threshold Model (based on Transformed wage)\n');
 [theta, vc, stderr]=tobit4(data,dep2,bin,ind,theta);
 
 %CI
 BCo=theta(2);
 BErr=stderr(2);
 LBCI=BCo-1.96*BErr;   %Lower CI
 HBCI=BCo+1.96*BErr;   %Upper CI
 fprintf('---Confidence Interval---\n');
 fprintf('(%2.4f,%2.4f)',LBCI, HBCI);
 
 %Part D
 fprintf('\n');
 fprintf('\n-----------Part D--------------- \n');
 
 %Backup Relevant Files
 theta2=theta;
 vc2=vc;
 
 %Create Geometric scaled wages
 load martins2;
 gmwage=geomean(wage2)
 gswage=(1/gmwage)*wage2;
 lgwage=log(gswage);
 save martins3;
 
 %define text strings for tobit models;
 data='martins3';
 dep2='gswage';
 ind='edu     pexp    pexp2   pexpchd pexpch2 ';
 bin='employ';
 
 fprintf('\n ---Geometrically Scaled Wages--- \n');
 
 % Use HL to generate startingvalues;
 [s1b,s2b]=tobit_hl2(data,dep2,bin,ind);
 theta=tobitPD(data,dep2,bin,ind,s2b);
 theta=tobitPD(data,dep2,bin,ind,theta);
 theta=tobitPD(data,dep2,bin,ind,theta);
 theta(7)=exp(theta(7));
 fprintf('\n Common Threshold Model (based on Transformed swage)\n');
 [theta, vc, stderr]=tobit4(data,dep2,bin,ind,theta);
 
 %define text strings for tobit models;
 data='martins3';
 dep2='lgwage';
 ind='edu     pexp    pexp2   pexpchd pexpch2 ';
 bin='employ';
 
 fprintf('\n ----------- \n');
 
 fprintf('\n ---Geometrically Scaled Log-Wages--- \n');
 
 % Use HL to generate startingvalues;
 [s1b,s2b]=tobit_hl2(data,dep2,bin,ind);
 theta=tobitPD(data,dep2,bin,ind,s2b);
 theta=tobitPD(data,dep2,bin,ind,theta);
 theta=tobitPD(data,dep2,bin,ind,theta);
 theta(7)=exp(theta(7));
 fprintf('\n Common Threshold Model (based on Transformed lgwage)\n');
 [theta, vc, stderr]=tobit4(data,dep2,bin,ind,theta);
 
 %closing output
 diary off
