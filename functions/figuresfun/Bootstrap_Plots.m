function [caux,InferenceMSW,figureorder] = Bootstrap_Plots(n,p,horizons,confidence,RForm, SVARinp, figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect)
%Implements bootstrap-type inference and produces plots comparing bootstrap inference and the delta-method. 
%   -Syntax:
%       [caux,InferenceMSW] = Bootstrap_Plots(n,p,horizons,confidence,RForm,SVARinp,figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect)
%   -Inputs:
%      n:               Number of variables in the VAR model 
%      p:               Number of lags in the VAR model                                                    (1 x 1)
%      horizons:        Number of horizons for the Impulse Response Functions (IRFs)                       (1 x 1)  
%                       (does not include the impact horizon 0)    
%      confidence:      Value for the standard and weak-IV robust confidence set                           (1 x 1)
%      RForm:           Structure containing the reduced form parameters
%      SVARinp:         Structure containing ydata, z, & n
%      figureorder:     Figure number                                                                      (1 x 1)
%      Plugin:          Structure containing standard plug-in inference
%      InferenceMSW:    InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%      time:            Time unit for the dataset (e.g. year, month, etc.)                                 (String)
%      columnnames:     Vector with the names for the endogenous variables, in the same order as ydata     (1 x n)
%      savdir:          Directory where the figures generated will be saved                                (String)
%      direct:          Directory where TestScript.m is located                                            (String) 
%      dataset_name:    The name of the dataset used for generating the figures (used in the output label) (String)
%      IRFselect:       Indices for the variables that the user wants separate IRF plots for
%      cumselect:       Indices for the variables that the user wants cumulative IRF plots for
%
%   -Output:
%       caux:
%       InferenceMSW:  Structure containing the MSW weak-iv robust confidence interval
%       figureorder:   Scalar that keeps track of the number of figures                                    (1 x 1)   
%                      generated


%% 1) Comparison of "standard" bootstrap inference and the delta-method

%check inputs
addpath('functions/AuxFunctions');

Bootstrap_Plots_Check(n,p,horizons,confidence, figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect)

%Non-cumulative graphs 
figureorder = figureorder + 1; 

figure(figureorder)

plots.order     = [1:n];

caux            = norminv(1-((1-confidence)/2),0,1);

for iplot = 1:n
    
    if n > ceil(sqrt(n)) * floor(sqrt(n))
            
        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));
    
    else
        
        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));
        
    end
    
    plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on
    
    g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,1),':b'); hold on
    
    dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));
    
    lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));
    
    h1    = plot(0:1:horizons,dmub,'--b'); hold on
    
    g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,1),':b'); hold on
    
    h2    = plot(0:1:horizons,lmub,'--b'); hold on
    
    clear dmub lmub
    
    h3 = plot([0 horizons],[0 0],'black'); hold off
    
    xlabel(time)
    
    title(columnnames(iplot));
    
    xlim([0 horizons]);
    
    if iplot == 1
        
        legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
            'CS^{Plug-in}')
        
        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','southeast')
     
    end
    
            
end
            
singletitle('Non-Cumulative IRFs (Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);
    

%Cumulative graphs 
figureorder = figureorder + 1; 

figure(figureorder)

plots.order     = [1:n];

caux            = norminv(1-((1-confidence)/2),0,1);

for iplot = 1:n
    
    if n > ceil(sqrt(n)) * floor(sqrt(n))
            
        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));
    
    else
        
        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));
        
    end
    
    
    plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on
    
    g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,2),':b'); hold on
    
    dmub  =  Plugin.IRFcum(iplot,:) + (caux*Plugin.IRFstderrorcum(iplot,:));
    
    lmub  =  Plugin.IRFcum(iplot,:) - (caux*Plugin.IRFstderrorcum(iplot,:));
    
    h1    = plot(0:1:horizons,dmub,'--b'); hold on
    
    g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,2),':b'); hold on
    
    h2    = plot(0:1:horizons,lmub,'--b'); hold on
    
    clear dmub lmub
    
    h3 = plot([0 horizons],[0 0],'black'); hold off
    
    xlabel(time)
    
    title(columnnames(iplot));
    
    xlim([0 horizons]);
    
    if iplot == 1
        
        legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
            'CS^{Plug-in}')
        
        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','southeast')
     
    end
    
            
end

singletitle('Cumulative IRFs (Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

%% 2) Save the output and plots in ./Output/Mat and ./Output/Figs
 
%Check if the Output File exists, and if not create one.
 
if exist(savdir,'dir')==0
    
    mkdir('savdir')
        
end

mat = strcat(savdir,'/Mat');

if exist(mat,'dir')==0
    
    mkdir(mat)
        
end

figs = strcat(savdir, '/Figs'); 

if exist(figs,'dir')==0
    
    mkdir(figs)
        
end

%Saving noncumulative plots
cd(mat);
 
output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
               'bootstrap_', num2str(100*confidence));
 
save(strcat('IRF_SVAR',output_label,'.mat'));
 
figure(figureorder-1)
 
cd(figs);
 
print(gcf,'-depsc2',strcat('IRF_SVAR',output_label,'.eps'));
 
%Saving cumulative plots

cd(mat);
 
output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
               'bootstrap_', num2str(100*confidence));
 
save(strcat('IRF_SVAR_CUM',output_label,'.mat'));
 
figure(figureorder)
 
cd(figs);
 
print(gcf,'-depsc2',strcat('IRF_SVAR_CUM',output_label,'.eps'));
 
cd(direct);

%% 3)Comparison of "standard" bootstrap inference and the delta-method (Selected IRF) 

if length(IRFselect) ~= 1
    
    if isempty(IRFselect) == 0

        figureorder = figureorder + 1; 

        figure(figureorder)

        plots.order     = 1:length(IRFselect);

        caux            = norminv(1-((1-confidence)/2),0,1);

        for i = 1:length(IRFselect) 

            iplot = IRFselect(i);

            if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))

                subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));

            else

                subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));

            end

            plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

            g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,1),':b'); hold on

            dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

            lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

            h1    = plot(0:1:horizons,dmub,'--b'); hold on

            g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,1),':b'); hold on

            h2    = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 horizons],[0 0],'black'); hold off

            xlabel(time)

            title(columnnames(iplot));

            xlim([0 horizons]);

            if i == 1

                legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
                    'CS^{Plug-in}')

                set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                legend boxoff

                legend('location','southeast')

            end


        end
            
        singletitle('Selected Non-Cumulative IRFs (Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

    end
    
end

%% 4) Generating separate bootstrap inference and delta method comparison for selected IRF and saving them to different folder

if isempty(IRFselect) == 0

    plots.order     = 1:length(IRFselect);

    caux            = norminv(1-((1-confidence)/2),0,1);

    for i = 1:length(IRFselect) 

        iplot = IRFselect(i);

        figureorder = figureorder + 1; 

        figure(figureorder);

        plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on

        g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,1),':b'); hold on

        dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));

        lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));

        h1    = plot(0:1:horizons,dmub,'--b'); hold on

        g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,1),':b'); hold on

        h2    = plot(0:1:horizons,lmub,'--b'); hold on

        clear dmub lmub

        h3 = plot([0 horizons],[0 0],'black'); hold off

        xlabel(time)

        title(columnnames(iplot));

        xlim([0 horizons]);

        legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
                'CS^{Plug-in}')

        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        legend boxoff

        legend('location','southeast')

        %Check if the Output File exists, and if not create one.

        MatBootIRFselect = strcat(savdir, '/Mat/MatBootIRFselect');

        if exist(MatBootIRFselect,'dir')==0

            mkdir(MatBootIRFselect)

        end

        FigsBootIRFselect = strcat(savdir, '/Figs/FigsBootIRFselect');

        if exist(FigsBootIRFselect,'dir')==0

            mkdir(FigsBootIRFselect)

        end

        cd(MatBootIRFselect);

        output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                   'bootstrap_', num2str(100*confidence), '_', num2str(iplot));

        save(strcat('IRF_SVAR',output_label,'.mat'),...
            'InferenceMSW','Plugin','RForm','SVARinp');

        figure(figureorder)

        cd(FigsBootIRFselect);

        print(gcf,'-depsc2',strcat('IRF_SVAR',output_label,'.eps'));

        cd(direct);

    end    

    clear plots output_label labelstrs dtype

else
    
end

%% 5) Comparison of "standard" bootstrap inference and the delta-method (Selected Cumulative) 

if length(cumselect) ~= 1

    if isempty(cumselect) == 0

        figureorder = figureorder + 1; 

        figure(figureorder)

        plots.order     = 1:length(cumselect);

        caux            = norminv(1-((1-confidence)/2),0,1);

        for i = 1:length(cumselect) 

            iplot = cumselect(i);

            if length(cumselect) > ceil(sqrt(length(cumselect))) * floor(sqrt(length(cumselect)))

                subplot(ceil(sqrt(length(cumselect))),ceil(sqrt(length(cumselect))),plots.order(1,i));

            else

                subplot(ceil(sqrt(length(cumselect))),floor(sqrt(length(cumselect))),plots.order(1,i));

            end

            plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on

            g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,2),':b'); hold on

            dmub  =  Plugin.IRFcum(iplot,:) + (caux*Plugin.IRFstderrorcum(iplot,:));

            lmub  =  Plugin.IRFcum(iplot,:) - (caux*Plugin.IRFstderrorcum(iplot,:));

            h1    = plot(0:1:horizons,dmub,'--b'); hold on

            g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,2),':b'); hold on

            h2    = plot(0:1:horizons,lmub,'--b'); hold on

            clear dmub lmub

            h3 = plot([0 5],[0 0],'black'); hold off

            xlabel(time)

            title(columnnames(iplot))

            xlim([0 horizons]);

            if i == 1

                legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
                    'CS^{Plug-in}')

                set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                legend boxoff

                legend('location','southeast')

            end


        end
            
        singletitle('Selected Cumulative IRFs (Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

    else

    end
    
else
    
end

%% 6) Generating separate bootstrap inference and delta method comparison for selected cumulative IRF and saving them to different folder
 
if isempty(cumselect) == 0
    plots.order     = 1:length(cumselect);

    caux            = norminv(1-((1-confidence)/2),0,1);

    for i = 1:length(cumselect) 

        iplot = cumselect(i);

        figureorder = figureorder + 1; 

        figure(figureorder);

        plot(0:1:horizons,Plugin.IRFcum(iplot,:),'b'); hold on

        g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2,2),':b'); hold on

        dmub  =  Plugin.IRFcum(iplot,:) + (caux*Plugin.IRFstderrorcum(iplot,:));

        lmub  =  Plugin.IRFcum(iplot,:) - (caux*Plugin.IRFstderrorcum(iplot,:));

        h1    = plot(0:1:horizons,dmub,'--b'); hold on

        g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1,2),':b'); hold on

        h2    = plot(0:1:horizons,lmub,'--b'); hold on

        clear dmub lmub

        h3 = plot([0 horizons],[0 0],'black'); hold off

        xlabel(time)

        title(strcat('Cumulative','{ }',columnnames(iplot)))

        xlim([0 horizons]);

        legend('SVAR-IV Estimator',strcat('Bootstrap CS^{Plug-in} (',num2str(100*confidence),'%)'),...
                'CS^{Plug-in}')

        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        legend boxoff

        legend('location','southeast')

        %Check if the Output File exists, and if not create one.

        MatBootIRFCUMselect = strcat(savdir, '/Mat/MatBootIRFCUMselect');

        if exist(MatBootIRFCUMselect,'dir')==0

            mkdir(MatBootIRFCUMselect)

        end

        FigsBootIRFCUMselect = strcat(savdir, '/Figs/FigsBootIRFCUMselect');

        if exist(FigsBootIRFCUMselect,'dir')==0

            mkdir(FigsBootIRFCUMselect)

        end

        cd(MatBootIRFCUMselect);

        output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                   'bootstrap_', num2str(100*confidence), '_', num2str(iplot));

        save(strcat('IRF_SVAR_CUM',output_label,'.mat'),...
            'InferenceMSW','Plugin','RForm','SVARinp');

        figure(figureorder)

        cd(FigsBootIRFCUMselect);

        print(gcf,'-depsc2',strcat('IRF_SVAR_CUM',output_label,'.eps'));

        cd(direct);

    end

    clear plots output_label labelstrs dtype
    
else
    
end

end

