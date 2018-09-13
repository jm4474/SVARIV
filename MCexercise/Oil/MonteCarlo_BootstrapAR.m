%% 1) READ ME

%Last update: August 15th, 2018. 

%This script implements a Monte-Carlo study to analyze the finite-sample
%coverage of the confidence interval suggested by Montiel-Olea, Stock,
%and Watson (2018). The Monte Carlo design is explained in our paper.
%(see MC design 1)
%Please change the auxparamMC value and the confidence level in Section 1
%to generate the appropriate figure.

%The script takes approximately 361 seconds in MacBook Pro @ 2.4
%3.4 GHz Intel Core i5 with 16 GB in Memory (OSX High Sierra)

%MATLAB Version: 9.4.0. (R2018a)
%MATLAB License Number: 650045
%Operating System: Mac OS X  Version: 10.13.6 Build: 17F77 
%Java Version: Java 1.7.0_75-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
%----------------------------------------------------------------------------------------------------
%MATLAB                                                Version 9.4         (R2018a)
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
           
application = 'Oil';      %Name of this empirical application. This name will be used for creating and accessing folders           

horizons    = 20;         %Number of horizons for IRF analysis

scale       = -20;        %Scale

NWlags      = 0;          %Newey-West lags

auxparamMC...
            = (5.5)^.5;   %Controls the size of the first-stage in the MC
 
MC.NB       = 1000;       %Number of samples from the asymptotic distribution
 
MC.norm     = 1;          %Norm: normalizing variable  
     
% Set auxparaMC = 18.95^.5 to get an implied first-stage of 10.17 
% Set auxparaMC = 5^.5     to get an implied first-stage of 3.87

confidence  = 0.95; 

IRFselect   = [2,3];
% By default, the program generates a single figure with the MC coverage for ALL variables
% in the VAR. However, IRFselect allows the user to generate an indepedent
% figure displaying the coverage of some specific variables of interest. 
% The program also saves the the coverage plots for the variables in IRFselect as separate .eps files 

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the above vector will select the variables 
% "Log(1/1-AMTR)","Log income", "Log Real GDP", "Unemployment Rate"

cumselect    = [1];
% cumselect allows the user to generate the coverge for cumulative IRFs of
% interest

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

MC.Gamma    = RForm.Gamma; %estimator of E[z_t eta_t]; n x 1;

MC.alphaaux = RForm.Gamma(1,1);

B1          = MC.Gamma./MC.alphaaux;

sigmae1     = abs(MC.alphaaux)./((MC.Gamma'*(MC.Sigma^(-1))*MC.Gamma)^.5); 

Baux        = null(((MC.Sigma^(-1/2))*B1)'); %orthonormal basis for the ortho-
%complement of Sigma^{-1/2} B_1 as described in the paper. 

MC.B        = [B1,(MC.Sigma^(1/2))*Baux]; clear Baux B1
%This is the matrix B that will be used to generate data from the SVAR

MC.n        = size(MC.AL,1);

MC.D        = diag([sigmae1^2,ones(1, MC.n-1)],0); clear sigmae1

%% 4) Set-up the parameters of the measurement error model for the external IV

Z         = RForm.externalIV;

MC.WHat   = RForm.WHat;

MC.sigma2Gamma1...
          = MC.WHat(((MC.n^2)*(RForm.p))+1,((MC.n^2)*(RForm.p))+1);

MC.T      = size(RForm.eta,2); %(change, if you want a larger sample size)
                               %(default is: size(RForm.eta,2))

MC.muZ    = mean(Z);

MC.varZ   = var(Z);

MC.alpha  = MC.alphaaux;

MC.impliedcorr...
          = MC.alpha./((MC.varZ^.5).*MC.D(1,1));

%MC.impliedfirststage ...
%         = MC.T*((MC.alpha*MC.B(1,1))^2)./MC.sigma2Gamma1;
     
%MC.impliedfirststagelim...
%         = MC.T*((MC.alpha*MC.B(1,1))^2)./(MC.Sigma(1,1)*((MC.muZ)^2+MC.varZ));     

MC.sigmav = (MC.varZ-(((auxparamMC*MC.alpha).^2)/MC.D(1,1)))^.5;   
% Note that auxparamMC cannot be too large, otherwise
% the sigmav becomes the square root of a negative number

clear Z;
    
MC.p      = RForm.p; 

MC.Y0     = RForm.Y0; 

MC.horizons = horizons;

MC.x      = scale; 

clear RForm hori scale

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

seed    = load('seed/seedMay12_MC.mat'); %The seed file is in the Seed folder

seed    = seed.seed;

rng(seed); clear seed

%% The loop for the Monte Carlo Study starts here:

%% 7) Generate a time series of length T+burnout from the DGP we have constructed

MCdraws = 1000;

burnout = 500; %Number of observations that will be discarded to eliminate
%the influence of initial conditions

coverageMCBoots...
        = zeros(MC.n,MC.horizons+1,MCdraws,2);
%Stores logical value 1 if the true IRF was covered by
%the Boostrap inference generated by each data for each draw

coverageMCMSW...
        = zeros(MC.n,MC.horizons+1,MCdraws,2);
%Stores logical value 1 if the true IRF was covered by
%the MSW C.I. generated by each data for each draw

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
            =fliplr(MC.Y0');              %Y_{t-1},Y_{t-2}, ... Y_{t-p}

    DATAMCaux ...
            = zeros(MC.n,MC.T+burnout);   %Initialize the data

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


    YMC      = DATAMCaux(:,(burnout-MC.p)+1:end);

    ZMC      = EIVMCaux(1,(burnout-MC.p)+1:end);

    MCdata.Y = YMC';
    
    MCdata.Z = ZMC';

    clearvars -except MC MCdata MCdraws mcdraw coverageMCMSW coverageMCdmethod burnout IRFMC FirstStageMC auxparamMC NWlags confidence coverageMCBoots application columnnames dataset_name cumselect IRFselect

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
    RFormMC.Gamma = RFormMC.eta*MCdata.Z(MC.p+1:end,1)/(size(RFormMC.eta,2)); %n times 1
    
    RFormMC.n     = MC.n;
    
    RFormMC.p     = MC.p;

    %b) Estimation of What

    [RFormMC.WHatall,RFormMC.WHat,RFormMC.V] = CovAhat_Sigmahat_Gamma(MC.p,RFormMC.X,MCdata.Z(MC.p+1:end,1),RFormMC.eta, NWlags);                

    cd ..

    cd ..

    %% 11) Some definitions for the next sections

    vechSigma   = RFormMC.V * RFormMC.Sigma(:);

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

    test_aux        = zeros(grid_size, ndraws, RFormMC.n, MC.horizons+1,2); % 5th dimension is for cumulative and non-cumulative

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

    % Collect the plug-in estimators of the IRF to analyze its finite-sample
    % distribution (cumulative and non cumulative)
    for i = 1:RFormMC.n

        IRFMC.IRFplugin(i,:,mcdraw, 1)     = PluginMC.IRF(i,:); % non cumulative

        IRFMC.IRFplugin(i,:,mcdraw, 2)     = PluginMC.IRFcum(i,:); % cumulative

        % Collect also the delta-method standard errors (not used). Cumulative
        % and non cumulative

        IRFMC.dmethodstderr(i,:,mcdraw, 1) = PluginMC.IRFstderror(i,:); % non cumulative

        IRFMC.dmethodstderr(i,:,mcdraw, 2) = PluginMC.IRFstderrorcum(i,:); % cumulative

    end

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

    reject       = squeeze(reject);
    %reject(n x horizons +1 x 2)

    not_reject   = 1 - reject;

    coverageMCBoots(:,:,mcdraw,:) = not_reject;
    %coverageMCBoots(n,horizons + 1, mcdraws,2) = reject;

    %MSW

    %Coverage of MSW (using the coefficients of the quadratic equation described in the paper)
        for i = 1:MC.n

            %non-cumulative
            clear aux;

            aux  = (MC.IRFZ(i,:,1).^2).*InferenceMSWMC.ahat(i,:) +... 
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

plots.order = 1:MC.n;
           
graphcount  = 1;

figure(graphcount);

graphcount  = graphcount + 1;

MC.impliedfirststage ...
            = round(mean(FirstStageMC),2);

addpath('functions/figuresfun');

%non-cumulative
for iplot = 1:MC.n
    
    subplot(3,1,iplot)
    
    plot(0:20,mean(coverageMCBoots(iplot,:,:,1),3),'o'); hold on
    
    plot(0:20,mean(coverageMCMSW(iplot,:,:,1),3),'rx'); hold on
    
    axis([0 MC.horizons .8 1]); 
    
    xlabel('Months after the shock');
    
    ylabel('MC Coverage');
    
    title(columnnames(iplot)); 
    
    if iplot == 2
        
        legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))
        
        legend('location','southeast')

    
    end

end

title_master = strcat('MC Coverage (',num2str(MCdraws),' MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);

%cumulative

figure(graphcount);

graphcount  = graphcount + 1;
 
for iplot   = 1:MC.n
    
    subplot(3,1,iplot)
    
    plot(0:20,mean(coverageMCBoots(iplot,:,:,2),3),'o'); hold on
    
    plot(0:20,mean(coverageMCMSW(iplot,:,:,2),3),'rx'); hold on
    
    axis([0 MC.horizons .8 1]); 
    
    xlabel('Months after the shock');
    
    ylabel('MC Coverage');
    
    title(columnnames(iplot)); 
    
    if iplot == 2
        
        legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))
        
        legend('location','southeast')
    
    end

end

title_master = strcat('Cumulative MC Coverage (',num2str(MCdraws),' MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);

%% 18) Save Coverage Plot

outdir = strcat('Output/',application,'/MC/Boots');

if exist(outdir,'dir')
    
    cd(outdir);
  
else
    
    mkdir(outdir);
    
    cd(outdir);

end

namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

if exist(namdir,'dir')
    
    cd(namdir);
  
else
    
    mkdir(namdir);
    
    cd(namdir);

end

output_label = strcat(dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_confidence=',num2str(confidence),'.eps');

save_count   = 1;

figure(graphcount-1)

print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_Cumulative_', output_label));

save_count   = save_count + 1;

figure(graphcount-2)

print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_', output_label));

save_count   = save_count + 1;

save(strcat('MC_Coverage_',dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_confidence=',num2str(confidence),'.mat'));

cd ..

cd ..

cd ..

cd ..

cd ..

%% 19) Plot Selected IRFs and Cumulatives

% This section of the script plots selected non-cumulative graphs for
% the selected IRF variables and selected cumulative graphs for the 
% selected cumselect variables. These plots are simply displayed, but
% not automatically saved. Section 17 generates the individual plots, which
% are saved.
if length(IRFselect) > 1
    
    figure(graphcount)

    graphcount  = graphcount + 1;

    plots.order = 1:length(IRFselect);

    for i = 1:length(IRFselect) 

        iplot   = IRFselect(i);

        if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))

            subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));

        else

            subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));

        end

        plot(0:20,mean(coverageMCBoots(iplot,:,:,1),3),'o'); hold on

        plot(0:20,mean(coverageMCMSW(iplot,:,:,1),3),'rx'); hold off

        axis([0 MC.horizons .8 1]); 

        xlabel('Months after the shock');

        ylabel('MC Coverage');

        title(columnnames(iplot));

        if i == 1

            legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))

            legend('location','southeast')


        end

    end

    clear title_master;

    title_master = strcat('Selected MC Coverage (',num2str(MCdraws),' MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

    singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);

end

%Cumulative plots

if length(cumselect) > 1
    
    figure(graphcount)

    graphcount  = graphcount + 1;

    plots.order = 1:length(cumselect);

    for i = 1:length(cumselect) 

        iplot   = cumselect(i);

        if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))

            subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));

        else

            subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));

        end

        plot(0:20,mean(coverageMCBoots(iplot,:,:,2),3),'o'); hold on

        plot(0:20,mean(coverageMCMSW(iplot,:,:,2),3),'rx'); hold off

        axis([0 MC.horizons .8 1]); 

        xlabel('Months after the shock');

        ylabel('MC Coverage');

        title(columnnames(iplot));
        
        if i == 1

            legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))

            legend('location','southeast')


        end

    end

    title_master = strcat('Cumulative MC Coverage (',num2str(MCdraws),' MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

    singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);
    
end


%% 20) Generating separate IRF (cumulative and non) plots

plots.order    = 1:length(IRFselect);

for i = 1:length(IRFselect) 

    iplot      = IRFselect(i);
    
    figure(graphcount);
    
    graphcount = graphcount + 1;
    
    plot(0:20,mean(coverageMCBoots(iplot,:,:,1),3),'o'); hold on
        
    plot(0:20,mean(coverageMCMSW(iplot,:,:,1),3),'rx'); hold off
    
    axis([0 MC.horizons .8 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat(columnnames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

	legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))

    legend('location','southeast')

end

plots.order    = 1:length(cumselect);

for i = 1:length(cumselect) 

    iplot      = cumselect(i);
    
    figure(graphcount);
    
    graphcount = graphcount + 1;
    
    plot(0:20,mean(coverageMCBoots(iplot,:,:,2),3),'o'); hold on
    
    plot(0:20,mean(coverageMCMSW(iplot,:,:,2),3),'rx'); hold off

    axis([0 MC.horizons .8 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat('Cumulative', {' '}, columnnames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

    legend('Bootstrap AR',strcat('MSW C.I (',num2str(100*confidence),'%)'))

	legend('location','southeast')

end

%% 21) Save separate selected IRF (cumulative and non) plots

outdir = strcat('Output/',application,'/MC/Boots');

if exist(outdir,'dir')
    
    cd(outdir);
  
else
    
    mkdir(outdir);
    
    cd(outdir);

end

output_label = strcat(dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_confidence=',num2str(confidence),'.eps');

namdir       = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

if exist(namdir,'dir')
    
    cd(namdir);
  
else
    
    mkdir(namdir);
    
    cd(namdir);

end

for i = 1:(length(cumselect)) 
    
    figure(graphcount - i)
    
    print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_Separate_Cumulative_', output_label));
    
    save_count = save_count + 1;

end

for i = length(cumselect)+1:(length(IRFselect)+length(cumselect))
    
    figure(graphcount - i)

    print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_Separate_', output_label, num2str(graphcount-i)));
    
    save_count = save_count + 1;
    
end

cd ..
    
cd ..

cd ..

cd ..

cd ..

%% 22) Box&Whisker Type plot to summarize the finite-sample distribution of IRF coefficients

figure(graphcount)

graphcount = graphcount + 1;

plots.order = 1:MC.n;

addpath('functions/figuresfun');

for iplot = 1:MC.n

    if MC.n > ceil(sqrt(MC.n)) * floor(sqrt(MC.n))

        subplot(ceil(sqrt(MC.n)),ceil(sqrt(MC.n)),plots.order(1,iplot));

    else

        subplot(ceil(sqrt(MC.n)),floor(sqrt(MC.n)),plots.order(1,iplot));

    end
    
    IRFplugin = squeeze(IRFMC.IRFplugin(iplot,:,:,1))';
    
    for hori  = 1:MC.horizons+1
       
        scatter(hori*ones(MCdraws,1),IRFplugin(:,hori),'b','.'); hold on
        
        scatter(hori,median(IRFplugin(:,hori)),'r','o'); hold on
        
        scatter(hori,quantile(IRFplugin(:,hori),.975),'r','^'); hold on
        
        scatter(hori,quantile(IRFplugin(:,hori),.025),'r','v'); hold on
        
        scatter(hori,mean(IRFplugin(:,hori)),'r','*'); hold on
  
    end
    
    clear IRFplugin;
    
    plot(MC.IRFZ(iplot,:,1),'r'); hold on % non-cumulative
    
    plot(MC.IRFChol(iplot,:,1),'--r'); hold on % non-cumulative

    hold off
    
    xlabel('Months after the shock');
    
    ylabel(strcat(columnnames(iplot),{' '},'estimators'));
    
end
    
title_master = strcat('IRF Plug-in Estimators',{' '},num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

title_master = title_master{1};

singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);

legend({'True IRF','Cholesky','Plug-in IRFs','Median','.975 quantile','.0275 quantile','Mean'},'Location','northeast')

%Cumulative 

figure(graphcount)

graphcount  = graphcount + 1;

plots.order = 1:MC.n;

for iplot   = 1:MC.n

    if MC.n > ceil(sqrt(MC.n)) * floor(sqrt(MC.n))

        subplot(ceil(sqrt(MC.n)),ceil(sqrt(MC.n)),plots.order(1,iplot));

    else

        subplot(ceil(sqrt(MC.n)),floor(sqrt(MC.n)),plots.order(1,iplot));

    end
    
    IRFplugin = squeeze(IRFMC.IRFplugin(iplot,:,:,2))';

    for hori  = 1:MC.horizons+1
       
        scatter(hori*ones(MCdraws,1),IRFplugin(:,hori),'b','.'); hold on
        
        scatter(hori,median(IRFplugin(:,hori)),'r','o'); hold on
        
        scatter(hori,quantile(IRFplugin(:,hori),.975),'r','^'); hold on
        
        scatter(hori,quantile(IRFplugin(:,hori),.025),'r','v'); hold on
        
        scatter(hori,mean(IRFplugin(:,hori)),'r','*'); hold on
  
    end
    
    clear IRFplugin;
    
    plot(MC.IRFZ(iplot,:,2),'r'); hold on % cumulative
    
    plot(MC.IRFChol(iplot,:,2),'--r'); hold on % cumulative

    hold off
        
    xlabel('Months after the shock');
        
    ylabel(strcat(columnnames(iplot),{' '},'estimators'));
    
end

title_master = strcat('Cumulative IRF Plug-in Estimators',{' '},num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

title_master = title_master{1};

singletitle(title_master,'fontsize',16,'xoff',0,'yoff',.03);

legend({'True IRF','Cholesky','Plug-in IRFs','Median','.975 quantile','.0275 quantile','Mean'},'Location','northeast')

%% 23) Save Box&Whisker Type plots

outdir = strcat('Output/',application,'/MC/Boots');

if exist(outdir,'dir')
    
    cd(outdir);
  
else
    
    mkdir(outdir);
    
    cd(outdir);

end

namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

if exist(namdir,'dir')
    
    cd(namdir);
  
else
    
    mkdir(namdir);
    
    cd(namdir);

end

output_label = strcat(dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_confidence=',num2str(confidence),'.eps');

figure(graphcount - 1)

print(gcf,'-depsc2',strcat(num2str(save_count),'Box&Whisker_Cumulative', output_label));

save_count = save_count + 1;

figure(graphcount - 2)

print(gcf,'-depsc2',strcat(num2str(save_count),'Box&Whisker', output_label)); % non-cumulative is not marked in the title

save_count = save_count + 1;

cd ..

cd ..

cd ..

cd ..

cd ..
