%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This is a test version to see whether the SVARIV_Luigi function is working: Auguest 13th, 2018
% Comment: We have tested this script on a Macbook Pro 
 
direct = pwd;

%% 1) Set number of VAR Lags, Newey West lags and confidence level.
 
p           = 2; %Number of lags in the VAR model
 
confidence  = .68; %Value for the standard and weak-IV robust confidence set

% Define the variables in the SVAR
columnnames = [{'Log(1/1-AMTR)'}, ...
               {'Log Income'}, ...
               {'Log Real GDP'},...
               {'Unemployment Rate'},...
               {'Inflation'}, {'FFR'}, {'Log GOV'}, {'Log RSTPrices'}, {'DLOG RDEBT'}];

time        = 'Year';  % Time unit for the dataset (e.g. year, month, etc.

NWlags      = 8;  % Newey-West lags

norm        = 1; % Variable used for normalization

scale       = -1; % Scale of the shock

horizons    = 6; %Number of horizons for the Impulse Response Functions(IRFs)

IRFselect   = [1, 3, 5, 8];
% The program generates a figure with the IRFs for all variables. However,
% you can specify with the IRFselect (variables of interest) which variables
% you also want to generate indepedent figures. 

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the above vector will select the variables "AMTR,
% Log Real GDP, Inflation and DLOG RDEBT."

 
%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

cd(strcat(pwd,'/Data'));
 
    years                               = xlsread('Time',... 
                                          'A2:A60');
                                        % time variable. Here as a 59x1
                                        % vector
 
    ydata                               = xlsread('Data',...
                                          'A2:I60');  
                                        % all of the endogenous variables
                                        % used in the VAR
                                      
    z                                   = xlsread('ExtIV',...
                                          'A2:A60');
                                        %External IV. 59x1 vector.
                                      
cd ..
 
%% 3 Test

savdir = strcat(direct,'/Output');  %selected directory where the output files will be saved
 
addpath(strcat(direct,'/functions'));

[Plugin, InferenceMSW, Chol, RForm] = SVARIV_General(p, confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, time);
 
% A more in depth description of the function can be found within the
% function file. For clarity purposes, we briefly describe the function
% below:
%
% -Inputs:
%       p:           Number of lags in the VAR model                                                (1 times 1)                                          
%       confidence:  Value for the standard and weak-IV robust confidence set                       (1 times 1) 
%       ydata:       Endogenous variables from the VAR model                                        (T times n) 
%       z:           External instrumental variable                                                 (T times 1)
%       NWlags:      Newey-West lags                                                                (1 times 1)
%       norm:        Variable used for normalization                                                (1 times 1)
%       scale:       Scale of the shock                                                             (1 times 1)
%       horizons:    Number of horizons for the Impulse Response Functions(IRFs)                    (1 times 1)
%       savdir:      Directory where the figures generated will be saved                            (String)
%       columnnames: Vector with the names for the endogenous variables, in the same order as ydata (1 times n)
%       time:        Time unit for the dataset (e.g. year, month, etc.)                             (String)
%
% -Output:
%       PLugin: Structure containing standard plug-in inference
%       InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%       Chol: Cholesky IRFs



%% 4) "Standard" bootstrap-type inference based on samples from the asy dist.

seed            = load(strcat(direct,'/seed/seedMay12.mat')); 
    
seed            = seed.seed;

disp('-')

disp('7) This section samples from the asy dist of the reduced-form parameters to conduct "standard" inference.')

disp('Standard inference based on sampling from the asy. dist. takes only:')  

tic;

n            = RForm.n;  %number of variables in the VAR

T = size(ydata, 2); % Number of observations/time periods.
 
[~,InferenceMSW.bootsIRFs] = ...
                  Gasydistboots(seed, 1000, n, p, norm, scale, horizons, confidence, T,...
                  RForm.AL(:),RForm.V*RForm.Sigma(:),RForm.Gamma(:),...
                  RForm.WHatall,@IRFSVARIV);              

toc;

%% 5) Comparison of "standard" bootstrap inference and the delta-method

disp('-')

disp('8) This section compares inference based on sampling from the asy-dist with delta-method inference')

figure(3)

plots.order     = [1:n];

caux            = norminv(1-((1-confidence)/2),0,1);

for iplot = 1:n
    
    if n > ceil(sqrt(n)) * floor(sqrt(n))
            
        subplot(ceil(sqrt(n)),ceil(sqrt(n)),plots.order(1,iplot));
    
    else
        
        subplot(ceil(sqrt(n)),floor(sqrt(n)),plots.order(1,iplot));
        
    end
    
    
    plot(0:1:horizons,Plugin.IRF(iplot,:),'b'); hold on
    
    %[~,~] = jbfill(0:1:hori,InferenceMSW.bootsIRFs(iplot,:,2),...
        %InferenceMSW.bootsIRFs(iplot,:,1),[204/255 204/255 204/255],...
        %[204/255 204/255 204/255],0,0.5); hold on
    
    g1    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,2),':b'); hold on
    
    dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));
    
    lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));
    
    h1    = plot(0:1:horizons,dmub,'--b'); hold on
    
    g2    =  plot(0:1:horizons,InferenceMSW.bootsIRFs(iplot,:,1),':b'); hold on
    
    h2    = plot(0:1:horizons,lmub,'--b'); hold on
    
    clear dmub lmub
    
    h3 = plot([0 5],[0 0],'black'); hold off
    
    xlabel(time)
    
    title(columnnames(iplot));
    
    if iplot == 1
        
        legend('SVAR-IV Estimator',strcat('AsyDist Std C.I (',num2str(100*confidence),'%)'),...
            'D-Method C.I.')
        
        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','southeast')
     
    end
    
    xlim([0 horizons-1]);
            
end
