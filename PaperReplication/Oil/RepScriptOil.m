%% 0) Introduction

% This Script generates Figures 1A & 1B for the empirical 
% illustration discussed in the paper 
% "Inference in Structural Vector Autoregressions With an External instrument",
% by José L. Montiel Olea, James H. Stock, and Mark W. Watson
% Working Paper, Columbia University (2018).

% This version: September 20th, 2018
% Comment:      We have tested this function on 
%               1) An iMac @3.4 GHz Intel Core i5 (16 GB 2400 MHz DDR4)
%                  Running Matlab R2018a.
%                  
%               2) A MacBook pro @2.4 GHz Intel Core i7 (8 GB 1600 MHz DDR3)
%                  Running Matlab R2016b.
%                  This script runs in about 60 seconds.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%MAKE SURE TO CHANGE THE CONFIDENCE LEVEL TO 95% TO GENERATE FIGURE 1B.
%The script generates Figures 1A & 1B one at a time so make sure to run 
%the script again with confidence = 0.95.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic 

clear 

cd ..

cd ..

direct = pwd;
%% 1) Set number of VAR Lags, Newey West lags and confidence level.

fprintf('Section 1 sets values for: \n')

fprintf('-VAR Lags \n')
fprintf('-Newey West lags (if desired) \n')
fprintf('-Confidence level \n')
fprintf('-VAR variable names \n')
fprintf('-Time unit \n')
fprintf('-Variable used for normalization \n')
fprintf('-Scale of the shock \n')
fprintf('-Number of horizons for the IRFs. \n')

application = 'Oil';  % Name of this empirical application. This name will be used for creating and accessing folders

p           = 24;     %Number of lags in the VAR model
 
confidence  = .68;    %Confidence Level for the standard and weak-IV robust confidence set,
                        %This confidence level generates Figure 1A, change to 0.95 to generate Figure 1B!

% Define the variables in the SVAR
columnnames = [{'Percent Change in Global Crude Oil Production'}, ...
               {'Index of real economic activity'}, ...
               {'Real Price of Oil'}];

time        = 'Month';  % Time unit for the dataset (e.g. year, month, etc).

NWlags      = 0;  % Newey-West lags(if it is neccessary to account for time series autocorrelation)
                  % (set it to 0 to compute heteroskedasticity robust std errors)

norm        = 1; % Variable used for normalization

scale       = 1; % Scale of the shock

horizons    = 20; %Number of horizons for the Impulse Response Functions(IRFs)
                 %(does not include the impact or horizon 0)
                 
%% 2) Load data (saved in structure "data")
 %  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

disp('-')

disp('Section 2 loads the data and specifies the name of the dataset.')

cd(strcat(direct,'/Data/',application));

    ydata = xlsread('Data'); 
    %The frequency of this data is 1973:2 - 2007:12
    %The file data.txt was obtained directly from the AER website


    z    = xlsread('ExtIV');
    %The frequency of this data is 1973:2 - 2004:09
    %The .xls file was created by Jim Stock and Mark Watson and it 
    %contains a monthly measure of Kilian's [2008] instrument for the
    %oil supply shock. 
    
    years = xlsread('time');
    
dataset_name = 'OilData'; %The name of the dataset used for generating the figures (used in the output label)

%% 3) Definitions for the next sections 

disp('-')

disp('Sections 3-6 generate the IRFs and confidence intervals using the ')
disp('methodology described in Montiel Olea, Stock and Watson (2018)') 

cd(direct)


SVARinp.ydata = ydata;

SVARinp.Z = z;

SVARinp.n        = size(ydata,2); %number of columns(variables)

RForm.p          = p; %RForm_user.p is the number of lags in the model

%% 4) Least-squares, reduced-form estimation

addpath('functions/RForm');

%a) Estimation of (AL, Sigma) and the reduced-form innovations
[RForm.mu, ...
     RForm.AL, ...
     RForm.Sigma,...
     RForm.eta,...
     RForm.X,...
     RForm.Y]        = RForm_VAR(SVARinp.ydata, p);

    %b) Estimation of Gammahat (n times 1)

    RForm.Gamma      = RForm.eta*SVARinp.Z(p+1:end,1)/(size(RForm.eta,2));   %sum(u*z)/T. Used for the computation of impulse response.
    %(We need to take the instrument starting at period (p+1), because
    %we there are no reduced-form errors for the first p entries of Y.)

    %c) Add initial conditions and the external IV to the RForm structure

    RForm.Y0         = SVARinp.ydata(1:p,:);

    RForm.externalIV = SVARinp.Z(p+1:end,1);

    RForm.n          = SVARinp.n;

%% 5) Estimation of the asymptotic variance of A,Gamma

% Definitions

n            = RForm.n; % Number of endogenous variables

T            = (size(ydata,1)); % Number of observations (time periods)

d            = ((n^2)*p)+(n);     %This is the size of (vec(A)',Gamma')'

dall         = d+ (n*(n+1))/2;    %This is the size of (vec(A)',vec(Sigma), Gamma')'

%a) Covariance matrix for vec(A,Gammahat). Used
%to conduct frequentist inference about the IRFs. 
 
[RForm.WHatall,RForm.WHat,RForm.V] = ...
    CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,NWlags);                
 
%NOTES:
%The matrix RForm.WHatall is the covariance matrix of 
% vec(Ahat)',vech(Sigmahat)',Gamma')'
 
%The matrix RForm.WHat is the covariance matrix of only
% vec(Ahat)',Gamma')' 
 
% The latter is all we need to conduct inference about the IRFs,
% but the former is needed to conduct inference about FEVDs.

%% 6) Compute standard and weak-IV robust confidence set suggested in MSW
 
 %Apply the MSW function
 
addpath(strcat(direct,'/functions/Inference'));
 
[InferenceMSW,Plugin,Chol] = MSWfunction(confidence,norm,scale,horizons,RForm,1);

%% 7) Plot Figures 1A or 1B 

disp('-')

disp('Section 7 in this script plots the IRFs with different confidence intervals ')
disp('in Figures 1A & 1B') 

addpath('functions/figuresfun'); 

figureorder = 1;

if confidence == 0.68

        caux            = norminv(1-((1-confidence)/2),0,1);

        figure(figureorder); 
        
        figureorder = figureorder + 1;

        subplot(3,1,1)

        iplot = 1; 

        plot(0:1:horizons,Chol(iplot,:,2),'r'); hold on

        f1 = plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on        

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWuboundcum(iplot,:),...
                InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRFcum(iplot,:) + (caux*Plugin.IRFstderrorcum(iplot,:));

            lmub  =  Plugin.IRFcum(iplot,:) - (caux*Plugin.IRFstderrorcum(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on        

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f1,'LineWidth',3); hold off

            xlabel(time)

            title(strcat('Cumulative',{' '},columnnames(iplot)));

            xlim([0 horizons]);

            axis([0 20 -2 2]);

            legend('Cholesky estimate','SVAR-IV estimate','CS^{AR}',...
                    'CS^{plug-in}')

            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            legend boxoff

            legend('location','southeast')

        subplot(3,1,2)

        iplot = 2; 

        f2 = plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
                InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

            lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:horizons,Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f2,'LineWidth',3); hold off

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

            axis([0 20 -0.5 2]);

        subplot(3,1,3)

        iplot = 3;

        f3 = plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
                InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

            lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:horizons,Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f3,'LineWidth',3); hold off

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

            axis([0 20 -1 3]);
            
            title = strcat('A. 68% Confidence Sets');
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);
            
            clear title;
            
            output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                       num2str(100*confidence));
                   
            figure(figureorder-1)

            cd('PaperReplication/Oil/Figures');

            print(gcf,'-depsc2',strcat('Figure1A',output_label,'.eps'));

            cd .. ;

elseif confidence == 0.95
        
        caux            = norminv(1-((1-confidence)/2),0,1);

        figure(figureorder); 
        
        figureorder = figureorder+1;

        subplot(3,1,1)

        iplot = 1; 

        plot(0:1:horizons,Chol(iplot,:,2),'r'); hold on

        f1 = plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on        

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWuboundcum(iplot,:),...
                InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRFcum(iplot,:) + (caux*Plugin.IRFstderrorcum(iplot,:));

            lmub  =  Plugin.IRFcum(iplot,:) - (caux*Plugin.IRFstderrorcum(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on        

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f1,'LineWidth',3); hold off

            xlabel(time)

            title(strcat('Cumulative',{' '},columnnames(iplot)));

            xlim([0 horizons]);

            axis([0 20 -2 2]); 

            legend('Cholesky estimate','SVAR-IV estimate','CS^{AR}',...
                    ' CS^{plug-in}')

            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            legend boxoff

            legend('location','southeast')

        subplot(3,1,2)

        iplot = 2; 

        f2 = plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
                InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

            lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:horizons,Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f2,'LineWidth',3); hold off

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

            axis([0 20 -0.5 2]); 

        subplot(3,1,3)

        iplot = 3;

        f3 = plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
                InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

            lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:horizons,dmub,'--b'); hold on

            h2 = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:horizons,Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 horizons],[0 0],'black'); hold off

            set(f3,'LineWidth',3); hold off

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

            axis([0 20 -1 3]); 
            
            title = strcat('B. 95% Confidence Sets');
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);
            
            clear title;
            
            output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                       num2str(100*confidence));
                   
            figure(figureorder-1)

            cd('PaperReplication/Oil/Figures');

            print(gcf,'-depsc2',strcat('Figure1B',output_label,'.eps'));

            cd .. ;
            
else
        
end

toc


%% 8) AR confidence set using bootstrap implementation

seed            = load(strcat(direct,'/seed/seedMay12.mat')); 
    
seed            = seed.seed;

NB = 1000;                      % Number of bootstrap replications 

ARylim(:,:,1) = [-1, 2; -0.5 1.5;-1 2];
ARylim(:,:,2) = [-1, 2; -0.5 1.5;-1 2];

cd(strcat(direct,'/functions/Inference'));

if confidence == 0.68

    multiplier = 1;         %Scalar that multiplies ARylim to incrase the lower and upper bounds for the grid of lambdas.

    grid_size = 1000;       % Number of lambdas in the grid, for each variable and for each horizon.

end

if confidence == 0.90

    multiplier = 1;         %Scalar that multiplies ARylim to incrase the lower and upper bounds for the grid of lambdas.

    grid_size = 1000;       % Number of lambdas in the grid, for each variable and for each horizon.

end

if confidence == 0.95

    multiplier = 1;        %Scalar that multiplies ARylim to incrase the lower and upper bounds for the grid of lambdas.

    grid_size = 1000;      % Number of lambdas in the grid, for each variable and for each horizon.

end

[reject, bootsIRFs, ~, ~, null_grid] = GasydistbootsAR(ydata, T, seed, RForm.n, NB, p, norm, scale, horizons, confidence, SVARinp, NWlags, RForm.AL, RForm.Sigma, RForm.Gamma, RForm.V, RForm.WHatall, Plugin, multiplier, grid_size,ARylim);
        
%% 9) AR Bootstrap plots

disp('-')

disp('Section 7 in this script plots the IRFs with different confidence intervals ')
disp('in Figures 1A & 1B') 

addpath('functions/figuresfun'); 

if confidence == 0.68

    caux            = norminv(1-((1-confidence)/2),0,1);

    figure(figureorder); 

    figureorder = figureorder + 1;

    subplot(3,1,1)

    iplot = 1; 

    plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on

    [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWuboundcum(iplot,:),...
        InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
        [204/255 204/255 204/255],0,0.5); hold on

    for hor = 0:horizons

        rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 1),iplot,2));

        unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 0),iplot,2));

        normalize = (iplot == norm && hor == 0);

        if (normalize == 1)

            l1 = plot(hor,1,'b--o');

        end

        if (isempty(rejected_grid) == 1 && normalize == 0)

            disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

        end

        if (isempty(unrejected_grid) == 0 && normalize == 0)

            %plot unrejected ones (also grouping points so that legend is
            %clean)
            l2 = plot(hor,unrejected_grid,'b--o'); hold on
            l2_group = hggroup; 
            set(l2,'Parent',l2_group)

        end

        clear rejected_grid unrejected_grid

    end

    h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

    xlabel(time)
    
    axis([0 20 -2 2]);

    title(strcat('Cumulative', {' '}, columnnames(iplot)));

    xlim([0 horizons]);
    
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
    set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

    legend('SVAR-IV Estimator','CS^{AR}',...
            'Bootstrap CS^{AR}')

    legend boxoff

    legend('location','southeast')
        
    for iplot = 2:3    
        
        subplot(3,1,iplot)
        
        plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

        [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
            InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on

        for hor = 0:horizons

            rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));

            unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1)); 

            normalize = (iplot == norm && hor == 0);

            if (normalize == 1)

                plot(hor,1,'b--o');

            end

            if (isempty(rejected_grid) == 1 && normalize == 0)

                disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot rejected ones
                %plot(hor,rejected_grid,'r--x'); hold on

                %plot %not rejected ones
                plot(hor,unrejected_grid,'b--o'); hold on

            end

            clear rejected_grid unrejected_grid
        
        end

        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

        xlabel(time)
        
        if iplot == 2
            
            axis([0 20 -0.5 2]);
            
        else
            
            axis([0 20 -1 3]);
            
        end

        title(columnnames(iplot));

        xlim([0 horizons]);
 
    end
    
    singletitle('Bootstrap CS^{AR} vs. CS^{AR}','fontsize',15,'xoff',0,'yoff',0.03);
    
    addpath('functions/figuresfun'); 

elseif confidence == 0.90

        caux            = norminv(1-((1-confidence)/2),0,1);

        figure(figureorder); 
        
        figureorder = figureorder + 1;

        subplot(3,1,1)

        iplot = 1; 
        
        plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on
    
        [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWuboundcum(iplot,:),...
            InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on




        for hor = 0:horizons
            
            rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 1),iplot,2));

            unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 0),iplot,2));
            
            normalize = (iplot == norm && hor == 0);

            if (normalize == 1)

                l1 = plot(hor,1,'b--o');

            end

            if (isempty(rejected_grid) == 1 && normalize == 0)

                disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot unrejected ones (also grouping points so that legend is
                %clean)
                l2 = plot(hor,unrejected_grid,'b--o'); hold on
                l2_group = hggroup; 
                set(l2,'Parent',l2_group)

            end

            clear rejected_grid unrejected_grid

        end
        
        axis([0 20 -2 2]);

        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

        xlabel(time)

        title(strcat('Cumulative', {' '}, columnnames(iplot)));

        xlim([0 horizons]);
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        
        legend('SVAR-IV Estimator','CS^{AR}',...
                'Bootstrap CS^{AR}')

        legend boxoff

        legend('location','southeast')
        
    for iplot = 2:3    
        
        subplot(3,1,iplot)
        
        plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

        [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
            InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on

        for hor = 0:horizons

            rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));

            unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1)); 

            normalize = (iplot == norm && hor == 0);

            if (normalize == 1)

                plot(hor,1,'b--o')

            end

            if (isempty(rejected_grid) == 1 && normalize == 0)

                disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot rejected ones
                %plot(hor,rejected_grid,'r--x'); hold on

                %plot %not rejected ones
                plot(hor,unrejected_grid,'b--o'); hold on

            end

            clear rejected_grid unrejected_grid

        end
        
        if iplot == 2
            
            axis([0 20 -0.5 2]);
            
        else
            
            axis([0 20 -1 3]);
            
        end

        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

        xlabel(time)

        title(columnnames(iplot));



        xlim([0 horizons]);
 
    end
    singletitle('Bootstrap CS^{AR} vs. CS^{AR}','fontsize',15,'xoff',0,'yoff',0.03);

elseif confidence == 0.95

        caux            = norminv(1-((1-confidence)/2),0,1);

        figure(figureorder); 
        
        figureorder = figureorder + 1;

        subplot(3,1,1)

        iplot = 1; 
        
        plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on
    
        [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWuboundcum(iplot,:),...
            InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on




        for hor = 0:horizons

            rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 1),iplot,2));
            
            unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,2) == 0),iplot,2));

            normalize = (iplot == norm && hor == 0);

            if (normalize == 1)

                l1 = plot(hor,1,'b--o');

            end

            if (isempty(rejected_grid) == 1 && normalize == 0)

                disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot unrejected ones (also grouping points so that legend is
                %clean)
                l2 = plot(hor,unrejected_grid,'b--o'); hold on
                l2_group = hggroup; 
                set(l2,'Parent',l2_group)

            end

            clear rejected_grid unrejected_grid

        end

        axis([0 20 -2 2]);

        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

        xlabel(time)

        title(strcat('Cumulative', {' '}, columnnames(iplot)));

        xlim([0 horizons]);
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

        legend('SVAR-IV Estimator','CS^{AR}',...
                'Bootstrap CS^{AR}')
            
        legend boxoff

        legend('location','southeast')
           
    
        
        for iplot = 2:3    
        
            subplot(3,1,iplot)

            plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
                InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on




            for hor = 0:horizons

                rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));

                unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1));            
            
                normalize = (iplot == norm && hor == 0);

                if (normalize == 1)

                    plot(hor,1,'b--o')

                end

                if (isempty(rejected_grid) == 1 && normalize == 0)

                    disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Try increasing the multiplier in GasydistbootsAR.m'));

                end

                if (isempty(unrejected_grid) == 0 && normalize == 0)

                    %plot rejected ones
                    %plot(hor,rejected_grid,'r--x'); hold on

                    %plot %not rejected ones
                    plot(hor,unrejected_grid,'b--o'); hold on

                end

                clear rejected_grid unrejected_grid

            end
            
            if iplot == 2

                axis([0 20 -0.5 2]);
                
            else

                axis([0 20 -1 3]);

            end
            
            h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

        end
        
        singletitle('Bootstrap CS^{AR} vs. CS^{AR}','fontsize',15,'xoff',0,'yoff',0.03);
    
end
   