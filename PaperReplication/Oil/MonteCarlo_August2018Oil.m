%% 1) READ ME

%Last update: August 15th, 2018. 

% This script replicates figures 2 from Montiel Olea, Stock and Watson (2018). It implements a Monte-Carlo study to analyze the finite-sample
% coverage of the confidence interval suggested by Montiel-Olea, Stock,
% and Watson (2016). The Monte Carlo design is explained in our paper.
%(see MC design 1)
% Please change the auxparamMC value in Section 1 to generate the
% appropriate figure.

% The script takes approximately 380 seconds in MacBook Pro @ 2.4
% 2.4 GHz Intel Core i7 with 8 GB in Memory (OSX High Sierra)

%MATLAB Version: 9.1.0.441655 (R2016b)
%MATLAB License Number: 650045
%Operating System: Mac OS X  Version: 10.13.5 Build: 17F77 
%Java Version: Java 1.7.0_75-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
%----------------------------------------------------------------------------------------------------
%MATLAB                                                Version 9.1         (R2016b)
%Simulink                                              Version 8.8         (R2016b)
%Bioinformatics Toolbox                                Version 4.7         (R2016b)
%Curve Fitting Toolbox                                 Version 3.5.4       (R2016b)
%Database Toolbox                                      Version 7.0         (R2016b)
%Econometrics Toolbox                                  Version 3.5         (R2016b)
%Financial Toolbox                                     Version 5.8         (R2016b)
%Global Optimization Toolbox                           Version 3.4.1       (R2016b)
%Mapping Toolbox                                       Version 4.4         (R2016b)
%Optimization Toolbox                                  Version 7.5         (R2016b)
%Statistics and Machine Learning Toolbox               Version 11.0        (R2016b)
%Symbolic Math Toolbox                                 Version 7.1         (R2016b)


cd ..

cd ..

addpath('functions/RForm');

% We suggest you to run this script section by section. 
beep off

clear; clc;

hori     = 20;   % Number of horizons for IRF analysis

x        = -20;  % Scale

NWlags   = 0;    % Newey-West lags

auxparamMC...
         = (5.5)^.5;   %Controls the size of the first-stage in the MC
 

% Set auxparaMC = 18.95^.5 to get an implied first-stage of 10.17 
% Set auxparaMC = 5^.5     to get an implied first-stage of 3.87


%% 2) Loading the main inputs for the Monte-Carlo design
%--------------------------------------
%(All the inputs are saved in the MCparameters folder)
%--------------------------------------  

RForm    = load('Data/Oil/RForm_Kilian_2018.mat');

MC.mu    = RForm.mu;  %vector of means; n x 1

MC.AL    = RForm.AL;  %vector of slope coefficients n x np

MC.Sigma = RForm.Sigma; %covariance matrix of residuals n x n


%% 3) Set-up the matrix B for the Monte-Carlo study

MC.Gamma = RForm.Gamma; %estimator of E[z_t eta_t]; n x 1;

MC.alphaaux...
         = (RForm.Gamma'*(RForm.Sigma^(-1))*RForm.Gamma)^(.5);
     
eaux     = [1;1;-1];    %Direction of B1 used in our design

B1       = eaux*((eaux'*(MC.Sigma)^(-1)*eaux)^(-.5)); %as described in our paper

clear eaux

Baux     = null(((MC.Sigma^(-1/2))*B1)'); %orthonormal basis for the ortho-
%complement of Sigma^{-1/2} B_1 as described in the paper. 

MC.B     = [B1,(MC.Sigma^(1/2))*Baux]; clear Baux B1
%This is the matrix B that will be used to generate data from the SVAR     

MC.D     = diag([1,1,1],0); 

    %% 4) Set-up the parameters of the measurement error model for the external IV

MC.n     = size(MC.AL,1);

Z        = RForm.externalIV;

MC.WHat   = RForm.WHat;

MC.sigma2Gamma1...
          = MC.WHat(((MC.n^2)*(RForm.p))+1,((MC.n^2)*(RForm.p))+1);

MC.T      = 356; %(change, if you want a larger sample size)
                               %(default is: size(RForm.eta,2))
MC.muZ  = mean(Z);

MC.varZ = var(Z);

MC.alpha = MC.alphaaux;

MC.impliedcorr...
         = MC.alpha./((MC.varZ^.5).*MC.D(1,1));

%MC.impliedfirststage ...
%         = MC.T*((MC.alpha*MC.B(1,1))^2)./MC.sigma2Gamma1;
     
%MC.impliedfirststagelim...
%         = MC.T*((MC.alpha*MC.B(1,1))^2)./(MC.Sigma(1,1)*((MC.muZ)^2+MC.varZ));     

MC.sigmav=(MC.varZ-(((auxparamMC*MC.alpha).^2)/MC.D(1,1)))^.5;     

clear Z;
    
MC.p    = RForm.p; 

MC.Y0   = RForm.Y0; 

MC.hori = hori;

MC.x    = x; 

clear RForm hori x

%% 5) Compute the true IRFs

%a) Compute the true MA coefficients

Caux    = [eye(MC.n),MARep(MC.AL,MC.p,MC.hori)]; 

C       = reshape(Caux,[MC.n,MC.n,MC.hori+1]);

Ccum    = cumsum(C,3);

B1      = MC.x*MC.B(:,1)./MC.B(1,1);

MC.IRFZ(:,:,1)...
        =reshape(sum(bsxfun(@times,C,B1'),2),[MC.n,MC.hori+1]);
    
MC.IRFZ(:,:,2)...
        =reshape(sum(bsxfun(@times,Ccum,B1'),2),[MC.n,MC.hori+1]);

%(:,:,1) Corresponds to noncumulative IRFs
%(:,:,2) Corresponds to cumulative IRFs

%In order to study the "bias" of the plug-in estimator we also report the
%Cholesky estimator

MC.B1Chol ...
        = MC.Sigma(:,1)./(MC.Sigma(1,1)^.5);
    
B1Chol  ...
        = MC.x*MC.B1Chol./MC.B1Chol(1,1);

MC.IRFChol(:,:,1)...
        =reshape(sum(bsxfun(@times,C,B1Chol'),2),[MC.n,MC.hori+1]);
    
MC.IRFChol(:,:,2)...
        =reshape(sum(bsxfun(@times,Ccum,B1Chol'),2),[MC.n,MC.hori+1]);

clear Caux C Ccum B1 B1chol;

%% 6) Set-up the residuals to generate a time series of the desired length from the SVAR-IV 

seed    = load('seed/seedMay12.mat'); %The seed file is in the Seed folder

seed    = seed.seed;

rng(seed); clear seed

%% The loop for the Monte Carlo Study starts here:

%% 7) Generate a time series of length T+burnout from the DGP we have constructed

MCdraws = 1000;

burnout = 500; %Number of observations that will be discarded to eliminate
%the influence of initial conditions

coverageMCMSW...
        = zeros(MC.n,MC.hori+1,MCdraws);
    
coverageMCdmethod...
        = zeros(MC.n,MC.hori+1,MCdraws);
    
IRFMC.IRFplugin...
        =zeros(MC.n,MC.hori+1,MCdraws);
    
IRFMC.dmethodstderr...
        =zeros(MC.n,MC.hori+1,MCdraws);

FirstStageMC...
        = zeros(1,MCdraws);

for mcdraw = 1:MCdraws   
    
lagsDATAMC...
        =zeros(MC.n,MC.p,1+MC.T+burnout);
    
lagsDATAMC(:,:,1)...
        =fliplr(MC.Y0'); %Y_{t-1},Y_{t-2}, ... Y_{t-p}
    
DATAMCaux ...
        = zeros(MC.n,MC.T+burnout);  %Initialize the data
    
EIVMCaux...
        = zeros(1,MC.T+burnout);      %Initialize the instrument
    
residuals...
        = randn(1+MC.n,MC.T+burnout); %Structural Innovations 
    
for ix=1:MC.T+burnout
    
    rformresiduals...
        = ((MC.B))*((MC.D)^.5)*residuals(1:MC.n,ix);
    
    DATAMCaux(1:MC.n,ix) ...
        = (MC.mu)+(MC.AL*reshape(lagsDATAMC(:,:,ix),[MC.n*MC.p,1]))...
        +(rformresiduals);
    
    lagsDATAMC(:,1,ix+1) ...
        = DATAMCaux(:,ix);
    
    lagsDATAMC(:,2:MC.p,ix+1) ...
        = lagsDATAMC(:,1:MC.p-1,ix);
    
    EIVMCaux(1,ix)...
        = MC.muZ+((auxparamMC*MC.alpha./(MC.D(1,1)^.5))*residuals(1,ix))+...
          (MC.sigmav*residuals(MC.n+1,ix));
    %Linear measurement error model for the external IV. Could be replaced
    %if desired.
    clear rformresiduals
end 

%% 8) Drop the first burnout-p observations (burn-out period)


YMC=DATAMCaux(:,(burnout-MC.p)+1:end);
ZMC=EIVMCaux(1,(burnout-MC.p)+1:end);

MCdata.Y=YMC';
MCdata.Z=ZMC';

clearvars -except MC MCdata MCdraws mcdraw coverageMCMSW coverageMCdmethod burnout IRFMC FirstStageMC auxparamMC NWlags

%Thus, an MC data set consists of two parts:
%i) The T times n vector YMC
%ii)The T times 1 vector ZMC

%These are the inputs required to construct the confidence interval
%in Montiel-Olea, Stock, and Watson 2016

%% 9) Use YMC in MC data to estimate the reduced-form parameters
%for the MC run.

cd('functions/RForm')

[RFormMC.mu,RFormMC.AL,RFormMC.Sigma,RFormMC.eta,RFormMC.X,RFormMC.Y] = RForm_VAR(MCdata.Y,MC.p);


%% 10) Use YMC and Z in MC data to estimate What

%a) Some definitions
RFormMC.Gamma= RFormMC.eta*MCdata.Z(MC.p+1:end,1)/(size(RFormMC.eta,2)); %n times 1
RFormMC.n=MC.n;
RFormMC.p=MC.p;

%b) Estimation of What

[RFormMC.WHatall,RFormMC.WHat,RFormMC.V] = CovAhat_Sigmahat_Gamma(MC.p,RFormMC.X,MCdata.Z(MC.p+1:end,1),RFormMC.eta, NWlags);                

cd ..

cd ..

%% 11) Use RForm.MC to estimate the MSW confidence interval
%(Note that analyzing coverage does not require computing the
%full confidence interval, but we do it to keep the code as
%simple as possible)

addpath('functions/Inference')

[InferenceMSWMC,PluginMC,~] = MSWfunction(.95,1,MC.x,MC.hori,RFormMC,0);

%Collect the plug-in estimators of the IRF to analyze its finite-sample
%distribution
IRFMC.IRFplugin(1,:,mcdraw)=PluginMC.IRFcum(1,:);
IRFMC.IRFplugin(2,:,mcdraw)=PluginMC.IRF(2,:);
IRFMC.IRFplugin(3,:,mcdraw)=PluginMC.IRF(3,:);

%Collect also the delta-method standard errors (not used)
IRFMC.dmethodstderr(1,:,mcdraw)=PluginMC.IRFstderrorcum(1,:);
IRFMC.dmethodstderr(2,:,mcdraw)=PluginMC.IRFstderror(2,:);
IRFMC.dmethodstderr(3,:,mcdraw)=PluginMC.IRFstderror(3,:);

%First-stage Stat

FirstStageMC(1,mcdraw) ...
      = (((InferenceMSWMC.T^.5)...
         *RFormMC.Gamma(1,1))^2)/...
         RFormMC.WHat(((RFormMC.n^2)*RFormMC.p)+1,((RFormMC.n^2)*RFormMC.p)+1);

clear MCdata RFormMC PluginMC

%% 12) Check if the true IRFs are covered (MSW and delta method)
%For the first period we take the cumulative IRFs

%a) We first compute the coverage for the cumulative responses of 
%world oil production to an oil supply shock. 

%Coverage of MSW (using the coefficients of the quadratic equation described in the paper)
aux1=(MC.IRFZ(1,:,2).^2).*InferenceMSWMC.ahatcum(1,:) +... 
      MC.IRFZ(1,:,2).*InferenceMSWMC.bhatcum(1,:) + InferenceMSWMC.chatcum(1,:);
aux1(1,1)=-inf; %The first entry is always covered. 
coverageMCMSW(1,:,mcdraw)=bsxfun(@le,aux1,0); %clear aux1

%Coverage of Dmethod
coverageMCdmethod(1,:,mcdraw)=(MC.IRFZ(1,:,2)<=InferenceMSWMC.Dmethoduboundcum(1,:))...
    .*(InferenceMSWMC.Dmethodlboundcum(1,:)<=MC.IRFZ(1,:,2));
coverageMCdmethod(1,1,mcdraw)=1;

%b) We now compute the coverage for the noncumulative response of
%global real activity 

%Coverage of MSW
aux2=(MC.IRFZ(2,:,1).^2).*InferenceMSWMC.ahat(2,:) +... 
      MC.IRFZ(2,:,1).*InferenceMSWMC.bhat(2,:) + InferenceMSWMC.chat(2,:);
coverageMCMSW(2,:,mcdraw)=(aux2<=0); clear aux2

%Coverage of Dmethod
coverageMCdmethod(2,:,mcdraw)=(MC.IRFZ(2,:,1)<=InferenceMSWMC.Dmethodubound(2,:))...
    .*(InferenceMSWMC.Dmethodlbound(2,:)<=MC.IRFZ(2,:,1));


%c) We now compute the coverage for the noncumulative response of
%the price of oil 

%Coverage MSW

aux3=(MC.IRFZ(3,:,1).^2).*InferenceMSWMC.ahat(3,:) +... 
      MC.IRFZ(3,:,1).*InferenceMSWMC.bhat(3,:) + InferenceMSWMC.chat(3,:);
coverageMCMSW(3,:,mcdraw)=(aux3<=0); clear aux3


%Coverage Dmethod
coverageMCdmethod(3,:,mcdraw)=(MC.IRFZ(3,:,1)<=InferenceMSWMC.Dmethodubound(3,:))...
    .*(InferenceMSWMC.Dmethodlbound(3,:)<=MC.IRFZ(3,:,1));

clearvars -except MC coverageMCMSW coverageMCdmethod MCdraws mcdraw InferenceMSWMC burnout IRFMC FirstStageMC auxparamMC NWlags

end
%% 12) Plot coverage 

MC.impliedfirststage ...
         = round(mean(FirstStageMC),2);

figure(1)
subplot(3,1,1)
plot(mean(coverageMCMSW(1,:,:),3),'o'); hold on
plot(mean(coverageMCdmethod(1,:,:),3),'rx'); hold off
axis([1 MC.hori .8 1]); 
xlabel('Months after the shock');
ylabel('MC Coverage');
title(strcat('Cumulative Response of Oil Production (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 
%%the command ¡°round(MC.alpha,2)¡± is not available for the MATLAB 2014a but
%%is viable and also recommended in later version, like MATLAB 2015a.
%%In older version, the user should use the ¡°roundn(X,-N)¡± syntax of the MATLAB
%%ROUND function instead. 
%%Note that the sign of N is reversed: ¡°round(X,N)¡± should be replaced with ¡°roundn(X,-N)¡±.


subplot(3,1,2)
plot(mean(coverageMCMSW(2,:,:),3),'o'); hold on
plot(mean(coverageMCdmethod(2,:,:),3),'rx'); hold off
axis([1 MC.hori .8 1]); 
xlabel('Months after the shock');
ylabel('MC Coverage');
legend('MSW','Delta-Method','Location','southeast');
title(strcat('Response of Global Real Activity (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 
 
subplot(3,1,3)
plot(mean(coverageMCMSW(3,:,:),3),'o'); hold on
plot(mean(coverageMCdmethod(3,:,:),3),'rx'); hold off
axis([1 MC.hori .8 1]); 
xlabel('Months after the shock');
ylabel('MC Coverage');
title(strcat('Response of the Real Price of Oil (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

%% 13) Box&Whisker Type plot to summarize the finite-sample distribution of IRF coefficients

figure(2)
subplot(3,1,1)
IRFplugin1=squeeze(IRFMC.IRFplugin(1,:,:))';
for hori=1:MC.hori+1
scatter(hori*ones(MCdraws,1),IRFplugin1(:,hori),'b','.'); hold on
scatter(hori,median(IRFplugin1(:,hori)),'r','o'); hold on
scatter(hori,quantile(IRFplugin1(:,hori),.975),'r','^'); hold on
scatter(hori,quantile(IRFplugin1(:,hori),.025),'r','v'); hold on
scatter(hori,mean(IRFplugin1(:,hori)),'r','*'); hold on
end
plot(MC.IRFZ(1,:,2),'r'); hold on
plot(MC.IRFChol(1,:,2),'--r'); hold on
hold off
axis([1 MC.hori+1 -40 40])
xlabel('Months after the shock');
ylabel('IRF Plug-in Estimators');
title(strcat('Cumulative Response of Oil Production (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

subplot(3,1,2)
IRFplugin2=squeeze(IRFMC.IRFplugin(2,:,:))';
plot(MC.IRFZ(2,:,1),'r'); hold on
plot(MC.IRFChol(2,:,1),'--r'); hold on
for hori=1:MC.hori+1
scatter(hori*ones(MCdraws,1),IRFplugin2(:,hori),'b','.'); hold on
scatter(hori,median(IRFplugin2(:,hori)),'r','o'); hold on
scatter(hori,quantile(IRFplugin2(:,hori),.975),'r','^'); hold on
scatter(hori,quantile(IRFplugin2(:,hori),.025),'r','v'); hold on
scatter(hori,mean(IRFplugin2(:,hori)),'r','*'); hold on
end
hold off
axis([1 MC.hori+1 -40 40])
legend('True IRF','Cholesky','Plug-in IRFs','Median','.975 quantile','.0275 quantile','Mean')
xlabel('Months after the shock');
ylabel('IRF Plug-in Estimators');
title(strcat('Response of Global Real Activity (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

subplot(3,1,3)
IRFplugin3=squeeze(IRFMC.IRFplugin(3,:,:))';
for hori=1:MC.hori+1
scatter(hori*ones(MCdraws,1),IRFplugin3(:,hori),'b','.'); hold on
scatter(hori,median(IRFplugin3(:,hori)),'r','o'); hold on
scatter(hori,quantile(IRFplugin3(:,hori),.975),'r','^'); hold on
scatter(hori,quantile(IRFplugin3(:,hori),.025),'r','v'); hold on
scatter(hori,mean(IRFplugin3(:,hori)),'r','*'); hold on
end
plot(MC.IRFZ(3,:,1),'r'); hold on
plot(MC.IRFChol(3,:,1),'--r'); hold on
hold off
axis([1 MC.hori+1 -40 40])
xlabel('Months after the shock');
ylabel('IRF Plug-in Estimators');
title(strcat('Response of the Real Price of Oil (',num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

%% 14) Save Coverage Plot

cd('PaperReplication/Oil/Figures/MC');
namdir=strcat(date,'_','K','_','p=',num2str(MC.p),'confidence0.95','_Design2');
if exist(namdir,'dir')
    cd(namdir);
    print(figure(1),'-depsc2',strcat('MC_Coverage_KilianDGP','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.eps'));
    print(figure(2),'-depsc2',strcat('MC_IRFDist','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.eps'));
    save(strcat('MC_Coverage_KilianDGP','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.mat'));
    cd ..
    cd ..
    cd ..
    cd ..
    cd ..
else    
    mkdir(namdir);
    cd(namdir);
    print(figure(1),'-depsc2',strcat('MC_Coverage_KilianDGP','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.eps'));
    print(figure(2),'-depsc2',strcat('MC_IRFDist','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.eps'));
    save(strcat('MC_Coverage_KilianDGP','_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.mat'));
    cd ..
    cd ..
    cd ..
    cd ..
    cd ..
end
