%% 0) Introduction

% This Script generates the figures for the empirical 
% illustration discussed in the paper 
% "Inference in Structural Vector Autoregressions with an external instrument",
% by José L. Montiel Olea, James H. Stock, and Mark W. Watson
% Working Paper, Columbia University (2018).

% This version: August 28th, 2018
% Comment:      We have tested this function on an iMac 
%               @3.4 GHz Intel Core i5 (16 GB 2400 MHz DDR4)
%               Running Matlab R2018a.
%               This script runs in about 1 second. 

tic 

%% 1) Load data for both 68% and 95% confidence intervals

% The data used in this script was generated using the SVARIV_General
% function and data from Kilian (2018). It
% contains the IRF estimates, confidence intervals and Cholesky estimates. 
% One can use the 'TestScriptOil.m' file to load Kilian's dataset and call
% SVARIV_General.

SForm_68 = load('Data/SForm_conf=0.68.mat');
SForm_95 = load('Data/SForm_conf=0.95.mat');

%% 2) Plots in new version of paper (run for 68% confidence sets)

% Section 2 plot the figures specific to the oil application (as used
% in the latest version of the paper). The paper plots one cumulative
% variable and two non-cumulative variables in one figure. The general
% version of this section (in SVARIV_General function section 3) plots all cumulative
% variables in one figure and all non-cumulative variables in another
% figure.

for confidence = [0.68,0.95]
    if confidence == 0.68

        caux            = norminv(1-((1-confidence)/2),0,1);

        figureorder = 1; 

        figure(figureorder); 

        subplot(3,1,1)

        iplot = 1; 

        plot(0:1:SForm_68.horizons,SForm_68.Chol(iplot,:,2),'r'); hold on

        f1 = plot(0:1:SForm_68.horizons,SForm_68.Plugin.IRFcum(iplot,:),'b'); hold on        

            [~,~] = jbfill(0:1:SForm_68.horizons,SForm_68.InferenceMSW.MSWuboundcum(iplot,:),...
                SForm_68.InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_68.Plugin.IRFcum(iplot,:) + (caux*SForm_68.Plugin.IRFstderrorcum(iplot,:));

            lmub  =  SForm_68.Plugin.IRFcum(iplot,:) - (caux*SForm_68.Plugin.IRFstderrorcum(iplot,:));

            h1 = plot(0:1:SForm_68.horizons,dmub,'--b'); hold on        

            h2 = plot(0:1:SForm_68.horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 SForm_68.horizons],[0 0],'black'); hold off

            set(f1,'LineWidth',3); hold off

            xlabel(SForm_68.time)

            title(strcat('Cumulative',{' '},SForm_68.columnnames(iplot)));

            xlim([0 SForm_68.horizons]);

            axis([0 20 0 1.5]);

            legend('Cholesky','SVAR-IV Estimator',strcat('MSW C.I (',num2str(100*confidence),'%)'),...
                    'D-Method C.I.')

            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            legend boxoff

            legend('location','southeast')

        subplot(3,1,2)

        iplot = 2; 

        f2 = plot(0:1:SForm_68.horizons,SForm_68.Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:SForm_68.horizons,SForm_68.InferenceMSW.MSWubound(iplot,:),...
                SForm_68.InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_68.Plugin.IRF(iplot,:) + (caux*SForm_68.Plugin.IRFstderror(iplot,:));

            lmub  =  SForm_68.Plugin.IRF(iplot,:) - (caux*SForm_68.Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:SForm_68.horizons,dmub,'--b'); hold on

            h2 = plot(0:1:SForm_68.horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:SForm_68.horizons,SForm_68.Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 SForm_68.horizons],[0 0],'black'); hold off

            set(f2,'LineWidth',3); hold off

            xlabel(SForm_68.time)

            title(SForm_68.columnnames(iplot));

            xlim([0 SForm_68.horizons]);

            axis([0 20 -0.2 0.2]);

        subplot(3,1,3)

        iplot = 3;

        f3 = plot(0:1:SForm_68.horizons,SForm_68.Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:SForm_68.horizons,SForm_68.InferenceMSW.MSWubound(iplot,:),...
                SForm_68.InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_68.Plugin.IRF(iplot,:) + (caux*SForm_68.Plugin.IRFstderror(iplot,:));

            lmub  =  SForm_68.Plugin.IRF(iplot,:) - (caux*SForm_68.Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:SForm_68.horizons,dmub,'--b'); hold on

            h2 = plot(0:1:SForm_68.horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:SForm_68.horizons,SForm_68.Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 SForm_68.horizons],[0 0],'black'); hold off

            set(f3,'LineWidth',3); hold off

            xlabel(SForm_68.time)

            title(SForm_68.columnnames(iplot));

            xlim([0 SForm_68.horizons]);

            axis([0 20 -0.5 0.3]);
            
            title = strcat('A. 68% Confidence Sets');
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);
            
            clear title;

    else

    end

    if confidence == 0.95

        caux            = norminv(1-((1-confidence)/2),0,1);

        figureorder = figureorder + 1; 

        figure(figureorder); 

        subplot(3,1,1)

        iplot = 1; 

        plot(0:1:SForm_95.horizons,SForm_95.Chol(iplot,:,2),'r'); hold on

        f1 = plot(0:1:SForm_95.horizons,SForm_95.Plugin.IRFcum(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:SForm_95.horizons,SForm_95.InferenceMSW.MSWuboundcum(iplot,:),...
                SForm_95.InferenceMSW.MSWlboundcum(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_95.Plugin.IRFcum(iplot,:) + (caux*SForm_95.Plugin.IRFstderrorcum(iplot,:));

            lmub  =  SForm_95.Plugin.IRFcum(iplot,:) - (caux*SForm_95.Plugin.IRFstderrorcum(iplot,:));

            h1 = plot(0:1:SForm_95.horizons,dmub,'--b'); hold on

            h2 = plot(0:1:SForm_95.horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 SForm_95.horizons],[0 0],'black'); hold off

            set(f1,'LineWidth',3); hold off

            xlabel(SForm_95.time)

            title(strcat('Cumulative',{' '},SForm_95.columnnames(iplot)));

            xlim([0 SForm_95.horizons]);

            axis([0 20 -1 2]);

            legend('Cholesky','SVAR-IV Estimator',strcat('MSW C.I (',num2str(100*confidence),'%)'),...
                    'D-Method C.I.')

            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            legend boxoff

            legend('location','southeast')

        subplot(3,1,2)

        iplot = 2; 

        f2 = plot(0:1:SForm_95.horizons,SForm_95.Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:SForm_95.horizons,SForm_95.InferenceMSW.MSWubound(iplot,:),...
                SForm_95.InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_95.Plugin.IRF(iplot,:) + (caux*SForm_95.Plugin.IRFstderror(iplot,:));

            lmub  =  SForm_95.Plugin.IRF(iplot,:) - (caux*SForm_95.Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:SForm_95.horizons,dmub,'--b'); hold on

            h2 = plot(0:1:SForm_95.horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:SForm_95.horizons,SForm_95.Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 SForm_95.horizons],[0 0],'black'); hold off

            set(f2,'LineWidth',3); hold off

            xlabel(SForm_95.time)

            title(SForm_95.columnnames(iplot));

            xlim([0 SForm_95.horizons]);

            axis([0 20 -0.5 1.5]);

        subplot(3,1,3)

        iplot = 3;

        f3 = plot(0:1:SForm_95.horizons,SForm_95.Plugin.IRF(iplot,:),'b'); hold on

            [~,~] = jbfill(0:1:SForm_95.horizons,SForm_95.InferenceMSW.MSWubound(iplot,:),...
                SForm_95.InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
                [204/255 204/255 204/255],0,0.5); hold on

            dmub  =  SForm_95.Plugin.IRF(iplot,:) + (caux*SForm_95.Plugin.IRFstderror(iplot,:));

            lmub  =  SForm_95.Plugin.IRF(iplot,:) - (caux*SForm_95.Plugin.IRFstderror(iplot,:));

            h1 = plot(0:1:SForm_95.horizons,dmub,'--b'); hold on

            h2 = plot(0:1:SForm_95.horizons,lmub,'--b'); hold on

            clear dmub lmub

            plot(0:1:SForm_95.horizons,SForm_95.Chol(iplot,:,1),'r'); hold on

            h3 = plot([0 SForm_95.horizons],[0 0],'black'); hold off

            set(f3,'LineWidth',3); hold off

            xlabel(SForm_95.time)

            title(SForm_95.columnnames(iplot));

            xlim([0 SForm_95.horizons]);

            axis([0 20 -1 2]);
            
            title = strcat('B. 95% Confidence Sets');
            
            singletitle(title,'fontsize',16,'xoff',0,'yoff',.03);


    else

    end
    
end

toc