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

%We suggest you to run this script section by section. 
beep off

clear; clc;

addpath('functions/RForm');

% Define the variables in the SVAR
columnnames = [{'Percent Change in Global Crude Oil Production'}, ...
               {'Index of real economic activity'}, ...
               {'Real Price of Oil'}];
           
application = 'Oil';       %Name of this empirical application. This name will be used for creating and accessing folders           

horizons   = 20;         %Number of horizons for IRF analysis

x          = -20;        %Scale

NWlags     = 0;          %Newey-West lags

auxparamMC...
           = (5.5)^.5;   %Controls the size of the first-stage in the MC
 
MC.NB      = 1000;       %Number of samples from the asymptotic distribution

MC.norm    = 1;          %Norm: normalizing variable  
     
% Set auxparaMC = 18.95^.5 to get an implied first-stage of 10.17 
% Set auxparaMC = 5^.5     to get an implied first-stage of 3.87

confidence = 0.95; 

dataset_name = 'OilData'; %The name of the dataset used for generating the figures (used in the output label)



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

MC.horizons = horizons;

MC.x    = x; 

clear RForm hori x

%% 5) Compute the true IRFs

%a) Compute the true MA coefficients

addpath('functions/StructuralIRF');

Caux    = [eye(MC.n),MARep(MC.AL,MC.p,MC.horizons)]; 

C       = reshape(Caux,[MC.n,MC.n,MC.horizons+1]);

Ccum    = cumsum(C,3);

B1      = MC.x*MC.B(:,1)./MC.B(1,1);

MC.IRFZ(:,:,1)...
        =reshape(sum(bsxfun(@times,C,B1'),2),[MC.n,MC.horizons+1]);
    
MC.IRFZ(:,:,2)...
        =reshape(sum(bsxfun(@times,Ccum,B1'),2),[MC.n,MC.horizons+1]);

%(:,:,1) Corresponds to noncumulative IRFs
%(:,:,2) Corresponds to cumulative IRFs

%In order to study the "bias" of the plug-in estimator we also report the
%Cholesky estimator

MC.B1Chol ...
        = MC.Sigma(:,1)./(MC.Sigma(1,1)^.5);
    
B1Chol  ...
        = MC.x*MC.B1Chol./MC.B1Chol(1,1);

MC.IRFChol(:,:,1)...
        =reshape(sum(bsxfun(@times,C,B1Chol'),2),[MC.n,MC.horizons+1]);
    
MC.IRFChol(:,:,2)...
        =reshape(sum(bsxfun(@times,Ccum,B1Chol'),2),[MC.n,MC.horizons+1]);

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

coverageMCBoots...
        = zeros(MC.n,MC.horizons+1,MCdraws,2);
    
coverageMCMSW...
        = zeros(MC.n,MC.horizons+1,MCdraws,2);
    
IRFMC.IRFplugin...
        =zeros(MC.n,MC.horizons+1,MCdraws);
    
IRFMC.dmethodstderr...
        =zeros(MC.n,MC.horizons+1,MCdraws);

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


YMC = DATAMCaux(:,(burnout-MC.p)+1:end);

ZMC = EIVMCaux(1,(burnout-MC.p)+1:end);

MCdata.Y=YMC';
MCdata.Z=ZMC';

clearvars -except MC MCdata MCdraws mcdraw coverageMCMSW coverageMCdmethod burnout IRFMC FirstStageMC auxparamMC NWlags confidence coverageMCBoots application columnnames dataset_name application

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


%% 11) Some definitions for the next sections

vechSigma = RFormMC.V * RFormMC.Sigma(:);

%% 12) Make sure that Whatall is symmetric and positive semidefinite

WHatall     = RFormMC.WHatall;

dall        = size(WHatall,1);

WHatall     = (WHatall + WHatall')/2;
    
[aux1,aux2] = eig(WHatall);
    
WHatall     = aux1*max(aux2,0)*aux1'; 

%% 13) Generate Draws for the Bootsrap
% Centered at (vec(AL)', Gamma')'

gvar    = [mvnrnd(zeros(MC.NB,dall),(WHatall)/MC.T)',...
                     zeros(dall,1)];
    
Draws   = bsxfun(@plus,gvar,...
          [RFormMC.AL(:);vechSigma;RFormMC.Gamma(:)]);

k       = size(RFormMC.Gamma,1)/RFormMC.n;        

%% 14) Evaluate the parameter of interest  
%(Note that analyzing coverage does not require computing the
%full confidence interval, so instead of testing for the entire grid
%of null hypotheses (lambdas) like we do in the GasydistbootsAR.m function
%here we simply test whether the true IRF is covered or not (the true IRF
%is our only lambda).

ndraws        = size(Draws,2);
     
pdSigma       = zeros(1,ndraws);

addpath('functions/StructuralIRF');

RFormIRFBoots = zeros(RFormMC.n, MC.horizons + 1,ndraws,2); %4th dimension corresponds to non-cumulative and cumulative values.

AlphaBoots    = zeros(1, ndraws);

for idraws    = 1:ndraws

    %i) Generate the draws for AL 
        
    AL        = reshape(Draws(1:(RFormMC.n^2)*RFormMC.p,idraws),[RFormMC.n,RFormMC.n*RFormMC.p]);

    %ii) Draws from Sigma

    vechSigma = Draws((RFormMC.n^2)*RFormMC.p+1:(RFormMC.n^2)*RFormMC.p+(RFormMC.n*(RFormMC.n+1)/2),idraws);

    Sigma     = tril(ones(RFormMC.n),0);
    
    Sigma(Sigma==1) = vechSigma';

    Sigma     = Sigma + tril(Sigma,-1)';

    %Check if the draws are positive definite

    if min(eig(Sigma))>0

        pdSigma(1,idraws) = 1;

    else

        pdSigma(1,idraws) = 0;

    end

    %iii) Draws from Gamma

    Gamma = reshape(Draws(((RFormMC.n^2)*RFormMC.p)+(RFormMC.n*(RFormMC.n+1)/2)+1:end,idraws),...
      [RFormMC.n,k]);  

    [RFormIRFBoots(:,:,idraws,:), AlphaBoots(idraws)] = IRFSVARIV(AL,Sigma,Gamma,MC.horizons,MC.x,MC.norm);
    
    RFormIRFBoots(:,:,idraws,:) = RFormIRFBoots(:,:,idraws,:).*AlphaBoots(idraws);
    % RFormIRFBoots(n,MC.horizons+1,idraws,2)
    
    clear AL vechSigma Gamma 
    
end

grid_size       = 1;

test_aux      = zeros(grid_size, ndraws, RFormMC.n, MC.horizons+1,2); % 5th dimension is for cumulative and non-cumulative

for var         = 1:RFormMC.n

    for horizon = 1:MC.horizons+1
        
        null_hyp = reshape(MC.IRFZ(var,horizon,:),[1,2]);
        % MC.IRFZ(3x21x2)
        
        %null_vec = reshape(null_grid(var,horizon,:),[grid_size,1]);
        
        test_aux(:,:,var,horizon,:) = ARTestStatistic(var, horizon, RFormIRFBoots, AlphaBoots, null_hyp, MC.T, ndraws);
        
    end
    
end

%recentering

AR_test         = (test_aux - test_aux(:,ndraws,:,:,:));
%grid_size, ndraws, RFormMC.n, MC.horizons+1,2

%% 15) Use RForm.MC to estimate the MSW confidence interval
%(Note that analyzing coverage does not require computing the
%full confidence interval, but we do it to keep the code as
%simple as possible)

addpath('functions/Inference')

[InferenceMSWMC,PluginMC,~] = MSWfunction(.95,1,MC.x,MC.horizons,RFormMC,0);

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

%% 16)  Check if the true IRFs are covered (Bootstrap Inference and MSW)

%Bootstrap Inference

aux          = reshape(pdSigma,[1,ndraws]);

bootsIRFs    = quantile(AR_test(:,aux==1,:,:,:),...
                          [((1-confidence)/2),1-((1-confidence)/2)],2);      
% bootsIRFs(grid_size, confidence intervals, variables, horizons+1, cumulative)

test_T = test_aux(:,ndraws,:,:,:);
% test_T(grid_size, 1, n, horizons+1,cumulative)

reject       = (test_T(:,:,:,:,:) < bootsIRFs(:,1,:,:,:)) | (test_T(:,:,:,:,:) > bootsIRFs(:,2,:,:,:));
% reject(grid_size, variables, horizons+1, cumulative)

reject = squeeze(reject);
%reject(n x horizons +1 x 2)

not_reject = 1 - reject;

coverageMCBoots(:,:,mcdraw,:) = not_reject;
%coverageMCBoots(n,horizons + 1, mcdraws,2) = reject;

%MSW

%Coverage of MSW (using the coefficients of the quadratic equation described in the paper)
    for i = 1:MC.n

        %non-cumulative
        clear aux;

        aux = (MC.IRFZ(i,:,1).^2).*InferenceMSWMC.ahat(i,:) +... 
              MC.IRFZ(i,:,1).*InferenceMSWMC.bhat(i,:) + InferenceMSWMC.chat(i,:);

        if i == 1

            aux(1,1) = -inf; %The first entry is always covered. 

        end    

        coverageMCMSW(i,:,mcdraw,1) = bsxfun(@le,aux,0); %clear aux1

        %cumulative
        clear aux_cum

        aux_cum = (MC.IRFZ(i,:,2).^2).*InferenceMSWMC.ahatcum(i,:) +... 
              MC.IRFZ(i,:,2).*InferenceMSWMC.bhatcum(i,:) + InferenceMSWMC.chatcum(i,:);

        if i == 1

            aux_cum(1,1) = -inf; %The first entry is always covered. 

        end    

        coverageMCMSW(i,:,mcdraw,2) = bsxfun(@le,aux_cum,0); %clear aux1

    end 

end

%% 17) Plot coverage 

plots.order     = 1:MC.n;
           
graphcount = 1;

figure(graphcount);

graphcount = graphcount + 1;

MC.impliedfirststage ...
         = round(mean(FirstStageMC),2);

addpath('functions/figuresfun');

%1) Cumulative response of Oil Production
    
subplot(3,1,1)

plot(0:20,mean(coverageMCBoots(1,:,:,2),3),'o'); hold on

plot(0:20,mean(coverageMCMSW(1,:,:,2),3),'rx'); hold on

axis([0 MC.horizons 0.8 1]); 

xlabel('Months after the shock');

ylabel('MC Coverage');

title(strcat('Cumulative', {' '},columnnames(1))); 


%2) Response of Global Real Activity
subplot(3,1,2)
    
plot(0:20,mean(coverageMCBoots(2,:,:,1),3),'o'); hold on

plot(0:20,mean(coverageMCMSW(2,:,:,1),3),'rx'); hold on

axis([0 MC.horizons 0.8 1]); 

xlabel('Months after the shock');

ylabel('MC Coverage');

title(columnnames(2)); 

legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))

legend('location','southeast')

%3) Response of the Real Price of Oil

subplot(3,1,3)

plot(0:20,mean(coverageMCBoots(3,:,:,1),3),'o'); hold on

plot(0:20,mean(coverageMCMSW(3,:,:,1),3),'rx'); hold on

axis([0 MC.horizons 0.5 0.8]);
 
xlabel('Months after the shock');

ylabel('MC Coverage');

title(columnnames(3)); 

title_master = strcat('MC Coverage (',num2str(MCdraws),' MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);


%% 18) Save Coverage Plot

outdir = strcat('Output/',application,'/MC/Boots');

if exist(outdir,'dir')
    
    cd(outdir);
  
else
    
    mkdir(outdir);
    
    cd(outdir);

end

output_label = strcat(dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_confidence=',num2str(confidence),'.eps');

save_count = 1;

namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

figure(graphcount-1)

print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_Cumulative_', output_label));

cd ..

cd ..

cd ..

cd ..




