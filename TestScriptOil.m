%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This is a test version to see whether the SVARIV_Luigi function is working: Auguest 13th, 2018
% Comment: We have tested this script on a Macbook Pro 
 
direct = pwd;

%% 1) Set number of VAR Lags, Newey West lags and confidence level.
 
fprintf('This script is an example on how to use the function SVARIV_General to report confidence intervals for  IRFs \n');

fprintf('estimated using the SVAR-IV approach described in MSW(18) and the function Gasydistboots.\n');

disp('-');

disp('Section 1 describes sets values for VAR Lags, Newey West lags and confidence level');

p           = 2; %Number of lags in the VAR model
 
confidence  = .68; %Confidence Level for the standard and weak-IV robust confidence set

% Define the variables in the SVAR
columnnames = [{'Percent Change in Global Crude Oil Production'}, ...
               {'Index of real economic activity'}, ...
               {'Real Price of Oil'}];

time        = 'Month';  % Time unit for the dataset (e.g. year, month, etc).

NWlags      = 8;  % Newey-West lags(if it is neccessary to account for time series autocorrelation)
                  % (set it to 0 to compute heteroskedasticity robust std errors)

norm        = 1; % Variable used for normalization

scale       = -1; % Scale of the shock

horizons    = 5; %Number of horizons for the Impulse Response Functions(IRFs)
                 %(does not include the impact or horizon 0)
                 
IRFselect   = [1, 2];
% By default, the program generates a single figure with the IRFs for ALL variables
% in the VAR. However, IRFselect allows the user to generate an indepedent
% figure displaying only some specific variables of interest. 
% The program also saves the IRFs in IRFselect as separate .eps files (To Do)

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the above vector will select the variables "AMTR,
% Log Real GDP, Inflation and DLOG RDEBT."

%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

disp('-')

disp('Section 2 in this script loads the data necessary for the SVAR-IV analysis.')

cd(strcat(pwd,'/5MSW_November_2016/2data'));

    ydata = xlsread('Data'); 
    %The frequency of this data is 1973:2 - 2007:12
    %The file data.txt was obtained directly from the AER website


    z    = xlsread('ExtIV');
    %The frequency of this data is 1973:2 - 2004:09
    %The .xls file was created by Jim Stock and Mark Watson and it 
    %contains a monthly measure of Kilian's [2008] instrument for the
    %oil supply shock. 
    
    years = xlsread('time');
    
    



cd(direct)
 
%% 3) Test

disp('-')

disp('Section 3 in this script calls the SVARIV_General function and samples from the asy dist of the reduced-form parameters to conduct "standard" inference.')

savdir = strcat(direct,'/Output');  %selected directory where the output files will be saved
 
addpath(strcat(direct,'/functions/Inference'));

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

disp('Section 4 in this script samples from the asy dist of the reduced-form parameters to conduct "standard" inference.')

disp('Standard inference based on sampling from the asy. dist. takes only:')  

tic;

n            = RForm.n;  %number of variables in the VAR

T  = size(ydata, 2); % Number of observations/time periods.

NB = 1000;           % Number of bootstrap repliactions 

[~,InferenceMSW.bootsIRFs] = ...
                  Gasydistboots(seed, 1000, n, p, norm, scale, horizons, confidence, T,...
                  RForm.AL(:),RForm.V*RForm.Sigma(:),RForm.Gamma(:),...
                  RForm.WHatall,@IRFSVARIV);              

toc;

%% 5) Comparison of "standard" bootstrap inference and the delta-method

disp('-')

disp('Finally, section 5 compares inference based on sampling from the asy-dist with delta-method inference')

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
    
    xlim([0 horizons]);
    
    if iplot == 1
        
        legend('SVAR-IV Estimator',strcat('AsyDist Std C.I (',num2str(100*confidence),'%)'),...
            'D-Method C.I.')
        
        set(get(get(g2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','southeast')
     
    end
    
            
end
