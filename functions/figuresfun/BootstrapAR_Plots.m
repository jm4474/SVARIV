function [] = BootstrapAR_Plots(n,p,horizons,confidence,RForm, SVARinp, figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect, reject, null_grid, norm)
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
%      reject:          4d logical array for whether an IRF is rejected or not, for                        (multiplier x n x horizons+1 x 2)
%                       each lambda, variable, horizon, cumulative and non-cumulative
%      null_grid:       grid of lambdas used in testing whether an IRF is                                  (grid_size x n x 2)
%                       rejected or not      
%      norm:            normalizing variable                                                               (1 x 1)
%
%   -Output:
%       caux:
%       InferenceMSW:  Structure containing the MSW weak-iv robust confidence interval

%% 1) Comparison of AR bootstrap inference and the MSW CI


%Non-cumulative graphs 
figureorder = figureorder + 1; 

figure(figureorder)

plots.order     = 1:n;

for iplot = 1:n

    if n > ceil(sqrt(n)) * floor(sqrt(n))
            
        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));
    
    else
        
        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));
        
    end
    
    plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on
    
    jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
        InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
        [204/255 204/255 204/255],0,0.5); hold on
        
  
    for hor = 0:horizons
        
        rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));
        
        unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1));
        
        normalize = (iplot == norm && hor == 0);
        
        if (normalize == 1)
                
            l1 = plot(hor,1,'b--o');
            
        end
            
        if (isempty(rejected_grid) == 1 && normalize == 0)
            
            disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));
        
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
    
    title(columnnames(iplot));

    xlim([0 horizons]);

    if iplot == 1
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        
        legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
            'Bootstrap CS-AR')

        legend boxoff
        
        legend('location','southeast')
     
    end

end
            
singletitle('Non-Cumulative IRFs (Weak IV Robust Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

%Cumulative graphs 
figureorder = figureorder + 1; 

figure(figureorder)

for iplot = 1:n
    
    if n > ceil(sqrt(n)) * floor(sqrt(n))
            
        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));
    
    else
        
        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));
        
    end
    
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
            
            disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));
        
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
    
    title(strcat(columnnames(iplot)));
    
    xlim([0 horizons]);

    if iplot == 1
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        
        legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
            'Bootstrap CS-AR')
        
        legend boxoff
        
        legend('location','southwest')
     
    end

end
            
singletitle('Cumulative IRFs (Weak IV Robust Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);
    
%% 2) Save the output and plots in ./Output/Mat and ./Output/Figs

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
               'ARbootstrap_', num2str(100*confidence));
 
save(strcat('IRF_SVAR',output_label,'.mat'));
 
figure(figureorder-1)
 
cd(figs);
 
print(gcf,'-depsc2',strcat('IRF_SVAR',output_label,'.eps'));
 
%Saving cumulative plots

cd(mat);
 
output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
               'ARbootstrap_', num2str(100*confidence));
 
save(strcat('IRF_SVAR_CUM',output_label,'.mat'));
 
figure(figureorder)
 
cd(figs);
 
print(gcf,'-depsc2',strcat('IRF_SVAR_CUM',output_label,'.eps'));
 
cd(direct);

%% 3)Comparison of AR bootstrap inference and the MSW CI (Selected IRF)

if length(IRFselect) ~= 1
    
    if isempty(IRFselect) == 0
        
        figureorder = figureorder + 1; 

        figure(figureorder)

        plots.order     = 1:length(IRFselect);

        for i = 1:length(IRFselect) 

            iplot = IRFselect(i);

            if length(IRFselect) > ceil(sqrt(length(IRFselect))) * floor(sqrt(length(IRFselect)))

                subplot(ceil(sqrt(length(IRFselect))),ceil(sqrt(length(IRFselect))),plots.order(1,i));

            else

                subplot(ceil(sqrt(length(IRFselect))),floor(sqrt(length(IRFselect))),plots.order(1,i));

            end
            
            plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on
            
            [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
            InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on
        
        
            for hor = 0:horizons
                
                rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));
        
                unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1));
                
                normalize = (iplot == norm && hor == 0);
        
                if (normalize == 1)
                
                    l1 = plot(hor,1,'b--o');
            
                end

                if (isempty(rejected_grid) == 1 && normalize == 0)

                    disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));

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

            title(columnnames(iplot));

            xlim([0 horizons]);

            if i == 1

                set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
                
                legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
                    'Bootstrap CS-AR')

                legend boxoff

                legend('location','southeast')

            end
            
        end
            
        singletitle('Selected Non-Cumulative IRFs (Weak IV Robust Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

    else
        
    end  
        
else

end

%% 4) Generating separate bootstrap inference and delta method comparison for selected IRF and saving them to different folder

if isempty(IRFselect) == 0

   
    for i = 1:length(IRFselect) 

        iplot = IRFselect(i);

        figureorder = figureorder + 1; 

        figure(figureorder);

        plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on
        
        [~,~] = jbfill(0:1:horizons,InferenceMSW.MSWubound(iplot,:),...
            InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
            [204/255 204/255 204/255],0,0.5); hold on


        for hor = 0:horizons
            
            rejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 1),iplot,1));
        
            unrejected_grid = squeeze(null_grid((reject(:,iplot,hor+1,1) == 0),iplot,1));


            normalize = (iplot == norm && hor == 0);

            if (normalize == 1)

                l1 = plot(hor,1,'b--o');

            end

            if (isempty(rejected_grid) == 1 && normalize == 0)

                disp(strcat('No values were rejected for the variable', {' '}, '"', columnnames(iplot),'"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot unrejected ones (also grouping points so that legend is
                %clean)
                 l2 = plot(hor,unrejected_grid,'b--o'); hold on
                 l2_group = hggroup; 
                 set(l2,'Parent',l2_group);

            end

            clear rejected_grid unrejected_grid

        end

        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

        xlabel(time)

        title(columnnames(iplot));

        xlim([0 horizons]);

        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

        legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
            'Bootstrap CS-AR')

        legend boxoff

        legend('location','southeast')
        
        %Check if the Output File exists, and if not create one.
        MatBootARIRFselect = strcat(savdir, '/Mat/MatARbootIRFselect');

        if exist(MatBootARIRFselect,'dir')==0

            mkdir(MatBootARIRFselect)

        end

        FigsBootARIRFselect = strcat(savdir, '/Figs/FigsARbootIRFselect');

        if exist(FigsBootARIRFselect,'dir')==0

            mkdir(FigsBootARIRFselect)

        end

        cd(MatBootARIRFselect);

        output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                   'ARbootstrap_', num2str(100*confidence), '_', num2str(iplot));

        save(strcat('IRF_SVAR',output_label,'.mat'),...
            'InferenceMSW','Plugin','RForm','SVARinp');

        figure(figureorder)

        cd(FigsBootARIRFselect);

        print(gcf,'-depsc2',strcat('IRF_SVAR',output_label,'.eps'));

        cd(direct);    
        
    end

else
    
end

%% 5) Comparison of AR bootstrap inference and the MSW CI (Selected Cumulative) 

if length(cumselect) ~= 1

    if isempty(cumselect) == 0

        figureorder = figureorder + 1; 

        figure(figureorder)

        plots.order     = 1:length(cumselect);

        for i = 1:length(cumselect) 

            iplot = cumselect(i);

            if length(cumselect) > ceil(sqrt(length(cumselect))) * floor(sqrt(length(cumselect)))

                subplot(ceil(sqrt(length(cumselect))),ceil(sqrt(length(cumselect))),plots.order(1,i));

            else

                subplot(ceil(sqrt(length(cumselect))),floor(sqrt(length(cumselect))),plots.order(1,i));

            end

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

                    disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));

                end

                if (isempty(unrejected_grid) == 0 && normalize == 0)
                    
                    %plot unrejected ones (also grouping points so that legend is
                    %clean)
                    l2 = plot(hor,unrejected_grid,'b--o'); hold on
                    l2_group = hggroup; 
                    set(l2,'Parent',l2_group);

                end
        
                clear rejected_grid unrejected_grid

            end
    
            h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0

            xlabel(time)

            title(strcat(columnnames(iplot)));

            xlim([0 horizons]);
            

            if i == 1

                set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

                legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
                    'Bootstrap CS-AR')

                legend boxoff

                legend('location','southeast')

            end
            
        end
            
        singletitle('Selected Cumulative IRFs (Weak IV Robust Bootstrap)','fontsize',16,'xoff',0,'yoff',0.04);

    else

    end
    
else
    
end

%% 6) Generating separate AR bootstrap inference and MSW CI comparison for selected cumulative IRF and saving them to different folder
 
if isempty(cumselect) == 0
    
    plots.order     = 1:length(cumselect);

    for i = 1:length(cumselect) 

        iplot = cumselect(i);

        figureorder = figureorder + 1; 

        figure(figureorder);

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

                disp(strcat('No values were rejected for the variable "Cumulative', {' '}, columnnames(iplot), '"', {' '}, 'horizon', {' '}, num2str(hor), '. Increase the multiplier in MSWfunction.m'));

            end

            if (isempty(unrejected_grid) == 0 && normalize == 0)

                %plot unrejected ones (also grouping points so that legend is
                %clean)
                l2 = plot(hor,unrejected_grid,'b--o'); hold on
                l2_group = hggroup; 
                set(l2,'Parent',l2_group);

            end
        
            clear rejected_grid unrejected_grid

        end
    
        h3 = plot([0 horizons],[0 0],'black'); hold off % black line at y = 0
    
        xlabel(time)
    
        title(strcat('Cumulative',{' '},columnnames(iplot)));
    
        xlim([0 horizons]);
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
        set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        set(get(get(l2_group,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');

        legend('SVAR-IV Estimator',strcat('CS-AR (',num2str(100*confidence),'%)'),...
            'Bootstrap CS-AR')

        legend boxoff

        legend('location','southeast')
        
        %Check if the Output File exists, and if not create one.

        MatBootARIRFCUMselect = strcat(savdir, '/Mat/MatARbootIRFCUMselect');

        if exist(MatBootARIRFCUMselect,'dir')==0

            mkdir(MatBootARIRFCUMselect)

        end

        FigsBootARIRFCUMselect = strcat(savdir, '/Figs/FigsARbootIRFCUMselect');

        if exist(FigsBootARIRFCUMselect,'dir')==0

            mkdir(FigsBootARIRFCUMselect)

        end

        cd(MatBootARIRFCUMselect);

        output_label = strcat('_p=',num2str(p),'_',dataset_name,'_',...
                   'bootstrapAR_', num2str(100*confidence), '_', num2str(iplot));

        save(strcat('IRF_SVAR_CUM',output_label,'.mat'),...
            'InferenceMSW','Plugin','RForm','SVARinp');

        figure(figureorder)

        cd(FigsBootARIRFCUMselect);

        print(gcf,'-depsc2',strcat('IRF_SVAR_CUM',output_label,'.eps'));

        cd(direct);

    end

    clear plots output_label labelstrs dtype
    
else
    
end
