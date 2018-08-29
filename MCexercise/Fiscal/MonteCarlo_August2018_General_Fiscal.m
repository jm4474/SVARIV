%% 1) READ ME

%Last update: August 20th, 2018. 

%UPDATE labels for the plots

% This script implements a Monte-Carlo study to analyze the finite-sample
% coverage of the confidence interval suggested by Montiel-Olea, Stock,
% and Watson (2018). 

% The design in this Monte-Carlo exercise in based on the parameters
% estimated from a Tax SVAR, but I can easily be adapted to another
% application by renaming the variables in Section 1 and replacing the 
% reduced-form estimators in Section 2. 

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

application = 'Fiscal';

% Define the variables in the SVAR
varNames = [{'Log(1/1-AMTR)'}, ...
               {'Log Income'}, ...
               {'Log Real GDP'},...
               {'Unemployment Rate'},...
               {'Inflation'}, {'FFR'}, {'Log GOV'}, {'Log RSTPrices'}, {'DLOG RDEBT'}];

hori         = 5;   % Number of horizons for IRF analysis

scale        = -1;  % Scale of the shock

auxparamMC   = 1;   % Controls the size of the first-stage in the MC
     
NWlags       = 0;   % Newey-West lags(if it is neccessary to account for time series autocorrelation)
                    % (set it to 0 to compute heteroskedasticity robust std errors)

confidence   = 0.95;% Confidence Level

norm         = 1;   %variable defining the normalization 

dataset_name = 'ALL(PS2003)'; 
                    %Input name of dataset, this will be used when creating names of data files and plots                                        

IRFselect    = [1,2,3,4];
% By default, the program generates a single figure with the MC coverage for ALL variables
% in the VAR. However, IRFselect allows the user to generate an indepedent
% figure displaying the coverage of some specific variables of interest. 
% The program also saves the the coverage plots for the variables in IRFselect as separate .eps files 

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the above vector will select the variables 
% "Log(1/1-AMTR)","Log income", "Log Real GDP", "Unemployment Rate"

cumselect = [5,6];
% cumselect allows the user to generate the coverge for cumulative IRFs of
% interest

% Make sure to match the indices above to the variables in your
% dataset. 

%% 2) Loading the main inputs for the Monte-Carlo design
%--------------------------------------
%(All the inputs are saved in the MCparameters folder)
%--------------------------------------  

RForm    = load(strcat('Data/',application,'/Tax_RForm.mat'));
                      %Reduced-form parameters for the MC
                      %These can be obtained by running TestScript and
                      %saving the RForm structure. 

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

Z               = RForm.externalIV;

MC.WHat         = RForm.WHat;

MC.sigma2Gamma1 = MC.WHat(((MC.n^2)*(RForm.p))+1,((MC.n^2)*(RForm.p))+1);

MC.T            = size(RForm.eta,2);
      
MC.muZ          = mean(Z);

MC.varZ         = var(Z);

MC.alpha        = MC.alphaaux;

MC.impliedcorr  = MC.alpha./((MC.varZ^.5).*MC.D(1,1));

MC.sigmav       = (MC.varZ-(((auxparamMC*MC.alpha).^2)/MC.D(1,1)))^.5;     
                % Note that auxparamMC cannot be too large, otherwise
                % the sigmav becomes the square root of a negative number

clear Z;
    
MC.p            = RForm.p; 

MC.Y0           = RForm.Y0; 

MC.hori         = hori;

MC.x            = scale; 

clear RForm hori x

%% 5) Compute the true IRFs

%a) Compute the true MA coefficients

addpath('functions/StructuralIRF');

Caux    = [eye(MC.n),MARep(MC.AL,MC.p,MC.hori)]; 

C       = reshape(Caux,[MC.n,MC.n,MC.hori+1]);

Ccum    = cumsum(C,3);

B1      = MC.x*MC.B(:,1)./MC.B(1,1);

MC.IRFZ(:,:,1)...
        =reshape(sum(bsxfun(@times,C,B1'),2),[MC.n,MC.hori+1]);
    
MC.IRFZ(:,:,2)...
        =reshape(sum(bsxfun(@times,Ccum,B1'),2),[MC.n,MC.hori+1]);

% (:,:,1) Corresponds to noncumulative IRFs
% (:,:,2) Corresponds to cumulative IRFs

% In order to study the "bias" of the plug-in estimator we also report the
% Cholesky estimator

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

seed    = load('seed/seedMay12_MC.mat'); %The seed file is in the Seed folder

seed    = seed.seed;

rng(seed); clear seed

%% The loop for the Monte Carlo Study starts here:

%% 7) Generate a time series of length T+burnout from the DGP we have constructed

MCdraws = 1000;

burnout = 500; % Number of observations that will be discarded to eliminate
%the influence of initial conditions

coverageMCMSW...
        = zeros(MC.n,MC.hori+1,MCdraws, 2);
    
coverageMCdmethod...
        = zeros(MC.n,MC.hori+1,MCdraws, 2);
    
IRFMC.IRFplugin...
        = zeros(MC.n,MC.hori+1,MCdraws, 2);
    
IRFMC.dmethodstderr...
        = zeros(MC.n,MC.hori+1,MCdraws, 2);

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


    YMC      = DATAMCaux(:,(burnout-MC.p)+1:end);
    
    ZMC      = EIVMCaux(1,(burnout-MC.p)+1:end);

    MCdata.Y = YMC';
    
    MCdata.Z = ZMC';

    clearvars -except MC MCdata MCdraws mcdraw coverageMCMSW coverageMCdmethod burnout IRFMC FirstStageMC auxparamMC confidence NWlags norm dataset_name varNames IRFselect cumselect application

    %Thus, an MC data set consists of two parts:
    %i) The T times n vector YMC
    %ii)The T times 1 vector ZMC

    %These are the inputs required to construct the confidence interval
    %in Montiel-Olea, Stock, and Watson 2018

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

    [RFormMC.WHatall,RFormMC.WHat,RFormMC.V] = CovAhat_Sigmahat_Gamma(MC.p,RFormMC.X,MCdata.Z(MC.p+1:end,1),RFormMC.eta,NWlags);                

    cd ..

    cd ..

    %% 11) Use RForm.MC to estimate the MSW confidence interval
    % (Note that analyzing coverage does not require computing the
    % full confidence interval, but we do it to keep the code as
    % simple as possible)
    % We calculate the confidence interval for both cumulative and
    % non comulative

    addpath('functions/Inference')

    [InferenceMSWMC,PluginMC,~]     = MSWfunction(confidence,norm,MC.x,MC.hori,RFormMC,0);

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
    
    n = RFormMC.n;
    
    clear MCdata RFormMC PluginMC

    %% 12) Check if the true IRFs are covered (MSW and delta method)

    % We now compute the coverage for the cumulative and non-cumulative responses of 
    % each variable to the shock of interest. 
    
    for i = 1:n
        % Coverage of MSW (using the coefficients of the quadratic equation described in the paper)

        % non-cumulative
        aux1                          = (MC.IRFZ(i,:,1).^2).*InferenceMSWMC.ahat(i,:) +... 
                                         MC.IRFZ(i,:,1).*InferenceMSWMC.bhat(i,:) + InferenceMSWMC.chat(i,:);

        if i == 1
            
            aux1(1,1)                     = -inf; %The first entry is always covered
        
        else
            
            % continue
        
        end
        
        coverageMCMSW(i,:,mcdraw, 1)  = bsxfun(@le,aux1,0); clear aux1 %non cumulative

        % cumulative
        aux1_cum                      = (MC.IRFZ(i,:,2).^2).*InferenceMSWMC.ahatcum(i,:) +... 
                                         MC.IRFZ(i,:,2).*InferenceMSWMC.bhatcum(i,:) + InferenceMSWMC.chatcum(i,:); %cumulative                               

        if i == 1
            
            aux1_cum(1,1)                 = -inf; %The first entry is always covered. 
        
        else
            
            %continue
        
        end

        coverageMCMSW(i,:,mcdraw, 2)  = bsxfun(@le,aux1_cum,0); clear aux1_cum %cumulative

        % Coverage of Dmethod

        %non-cumulative
        coverageMCdmethod(i,:,mcdraw,1) = (MC.IRFZ(i,:,1)<=InferenceMSWMC.Dmethodubound(i,:))...
                                        .*(InferenceMSWMC.Dmethodlbound(i,:)<=MC.IRFZ(i,:,1));


       % coverageMCdmethod(i,1,mcdraw,1) = 1;

        %cumulative
        coverageMCdmethod(i,:,mcdraw,2) = (MC.IRFZ(i,:,2)<=InferenceMSWMC.Dmethoduboundcum(i,:))...
                                        .*(InferenceMSWMC.Dmethodlboundcum(i,:)<=MC.IRFZ(i,:,2));

       % coverageMCdmethod(i,1,mcdraw,2) = 1;                            
    end

    clearvars -except MC coverageMCMSW coverageMCdmethod MCdraws mcdraw InferenceMSWMC burnout IRFMC FirstStageMC auxparamMC confidence NWlags norm n dataset_name varNames IRFselect cumselect application
end
%% 13) Plot coverage 

MC.impliedfirststage ...
         = round(mean(FirstStageMC),2);

plots.order     = 1:n;
           
graphcount = 1;

for i = 1:2
    
    figure(graphcount)
    
    graphcount = graphcount + 1;
    
    for iplot = 1:n

        if n > ceil(sqrt(n)) * floor(sqrt(n))

            subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));

        else

            subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));

        end

        plot(mean(coverageMCMSW(iplot,:,:,i),3),'o'); hold on
        
        plot(mean(coverageMCdmethod(iplot,:,:,i),3),'rx'); hold off
        
        axis([1 MC.hori 0 1]);
        
        xlabel('Months after the shock');
        
        ylabel(strcat(varNames(iplot), 'MC Coverage'));
        
        
    end
    
    if i == 1
        
            title = strcat('MC Coverage', num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');
            %For non-cumulative plots, we do not include this label in the
            %title
            
            addpath('functions/figuresfun');
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);
            
            clear title;

    else
        
            addpath('functions/figuresfun');
            
            title = strcat('Cumulative', {' '},'MC Coverage', {' '},varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');
            
            title = title{1};
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);
            
            clear title;
    end
    
end

%% 14) Save Coverage Plot

cd(strcat('Output/',application,'/MC'));

namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

output_label = strcat(dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.eps');

save_count = 1;

if exist(namdir,'dir')
    
    cd(namdir);
  
else
    
    mkdir(namdir);
    
    cd(namdir);

end

figure(graphcount-1)

print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_Cumulative_', output_label));

save_count = save_count + 1;

figure(graphcount-2)

print(gcf,'-depsc2',strcat(num2str(save_count),'MC_Coverage_', output_label));

save_count = save_count + 1;

save(strcat('MC_Coverage_',dataset_name,'_p=',num2str(MC.p),'_T=',num2str(InferenceMSWMC.T),'_MCFirstStage=',num2str(round(MC.impliedfirststage,2)),'.mat'));

cd ..

cd ..

cd ..

cd ..

    
%% 15) Plot Selected IRFs and Cumulatives

% This section of the script plots selected non-cumulative graphs for
% the selected IRF variables and selected cumulative graphs for the 
% selected cumselect variables 

figure(graphcount)

graphcount = graphcount + 1;

plots.order     = 1:length(IRFselect);

for i = 1:length(IRFselect) 

    iplot = IRFselect(i);
    
    if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))
            
        subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));
    
    else
        
        subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));
        
    end
    
    plot(mean(coverageMCMSW(iplot,:,:,1),3),'o'); hold on
        
    plot(mean(coverageMCdmethod(iplot,:,:,1),3),'rx'); hold off

    axis([1 MC.hori .5 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat(varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

end

figure(graphcount)

graphcount = graphcount + 1;

plots.order = 1:length(cumselect);

for i = 1:length(cumselect) 

    iplot = cumselect(i);
    
    if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))
            
        subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));
    
    else
        
        subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));
        
    end
    
    plot(mean(coverageMCMSW(iplot,:,:,2),3),'o'); hold on
        
    plot(mean(coverageMCdmethod(iplot,:,:,2),3),'rx'); hold off

    axis([1 MC.hori .5 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat('Cumulative', {' '}, varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

end


%% 16) Generating separate IRF (cumulative and non) plots

plots.order     = 1:length(IRFselect);

for i = 1:length(IRFselect) 

    iplot = IRFselect(i);
    
    
    figure(graphcount);
    
    graphcount = graphcount + 1;
    
    plot(mean(coverageMCMSW(iplot,:,:,1),3),'o'); hold on
        
    plot(mean(coverageMCdmethod(iplot,:,:,1),3),'rx'); hold off
    
    axis([1 MC.hori .5 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat(varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

end

plots.order = 1:length(cumselect);


for i = 1:length(cumselect) 

    iplot = cumselect(i);
    
    figure(graphcount);
    
    graphcount = graphcount + 1;
    
    plot(mean(coverageMCMSW(iplot,:,:,2),3),'o'); hold on
        
    plot(mean(coverageMCdmethod(iplot,:,:,2),3),'rx'); hold off

    axis([1 MC.hori .5 1]);

    xlabel('Months after the shock');

    ylabel('MC Coverage');

    title(strcat('Cumulative', {' '}, varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')'));

end

%% 17) Save separate selected IRF (cumulative and non) plots

cd(strcat('Output/',application,'/MC'));


namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

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

%% 18) Box&Whisker Type plot to summarize the finite-sample distribution of IRF coefficients

figure(graphcount)

graphcount = graphcount + 1;

plots.order = 1:n;

addpath('functions/figuresfun');

for iplot = 1:n

    if n > ceil(sqrt(n)) * floor(sqrt(n))

        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));

    else

        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));

    end
    
    IRFplugin = squeeze(IRFMC.IRFplugin(iplot,:,:,1))';

    for hori=1:MC.hori+1
       
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
    
    %axis([1 MC.hori+1 -40 40])  
    
    xlabel('Months after the shock');
    
    %ylabel('IRF Plug-in Estimators');
    
    ylabel(strcat(varNames(iplot),{' '},'estimators'));
    
    %title(strcat(varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

    
end

title = strcat('IRF Plug-in Estimators',{' '},num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

title = title{1};

singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);

legend({'True IRF','Cholesky','Plug-in IRFs','Median','.975 quantile','.0275 quantile','Mean'},'Location','northeast')

%Cumulative 

figure(graphcount)

graphcount = graphcount + 1;

plots.order = 1:n;

for iplot = 1:n

    if n > ceil(sqrt(n)) * floor(sqrt(n))

        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));

    else

        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));

    end
    
    IRFplugin = squeeze(IRFMC.IRFplugin(iplot,:,:,2))';

    for hori=1:MC.hori+1
       
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
    
    %axis([1 MC.hori+1 -40 40])  
    
    xlabel('Months after the shock');
    
    %ylabel('IRF Plug-in Estimators');
    
    ylabel(strcat(varNames(iplot),{' '},'estimators'));
    
    %title(strcat('Cumulative', {' '},varNames(iplot),num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')')); 

    
end

title = strcat('Cumulative IRF Plug-in Estimators',{' '},num2str(MCdraws),'MC draws, T=',num2str(InferenceMSWMC.T),', MC First Stage=',num2str(round(MC.impliedfirststage,2)),')');

title = title{1};

singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);

legend({'True IRF','Cholesky','Plug-in IRFs','Median','.975 quantile','.0275 quantile','Mean'},'Location','northeast')

%% 19) Save Box&Whisker Type plots

cd(strcat('Output/',application,'/MC'));

namdir = strcat(date,'_','K','_','p=',num2str(MC.p),num2str(confidence));

if exist(namdir,'dir')
    
    cd(namdir);
  
else
    
    mkdir(namdir);
    
    cd(namdir);

end

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
    