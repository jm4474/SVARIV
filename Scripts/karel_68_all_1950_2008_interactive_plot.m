%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This version: July 11th, 2017
% Comment: We have tested this script on a Macbook Pro 
%         @2.4 GHz Intel Core i7 (8 GB 1600 MHz DDR3)
%         Running Matlab R2016b.
%         This script runs in about 10 seconds.

clear; clc;

cd ..

main_d = pwd;

cd(main_d);

disp('This script reports confidence intervals for IRFs estimated using the SVAR-IV approach')

disp('(created by Karel Mertens and Jose Luis Montiel Olea)')

disp('-')

disp('This version: July 2017')

disp('-')

disp('(We would like to thank Qifan Han and Jianing Zhai for excellent research assistance)')

%% 1) Set number of VAR Lags, Newey West lags, the sub-dataset and confidence level.

disp('-')

disp('1) The first section defines the number of VAR Lags, Newey West lags, the sub-dataset and confidence level that will be used for local projection.')

prompt1 = 'Please input the number of VAR lag (our baseline is 2).';

p = input(prompt1); 

prompt2 = 'Please input the number of Newey-West lags (our baseline is 8).';

nw      = input(prompt2);

ok1     = 0; 

ok2     = 0;

while ok1 == 0
    
    strs1        = {'ALL(BR2011)', 'ALL(PS2003)','TOP 1%','TOP 5%',...
               'TOP 10%','TOP 5%-1%','TOP 10%-5%','BOTTOM 99%','BOTTOM 90%'};
           
    [dtype, ok1] = listdlg('PromptString','Select the dataset','SelectionMode','Single',...
                      'ListString', strs1);
                  
end

while ok2 == 0
    
    strs2        = {'68%','90%','95%'};
    
    [ctype, ok2] = listdlg('PromptString','Select the confidence level','SelectionMode','Single',...
                      'ListString', strs2);
                  
end

% Confidence level

if ctype == 1
    
        confidence = .68;
    
elseif ctype == 2
    
        confidence = .90;
        
else
    
        confidence = .95;
        
end

%% 2) Load data (saved in structure "data")

disp('-')

disp('2) The second section loads the data and saves it in the "data" structure')

cd(strcat(main_d,'/Data'));

[data.years, ~, ~]                      = xlsread('DATA_Mertens2015',... 
                                          'AMTR (Figure 1)', 'A6:A64');
                                      
if dtype == 1
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'B6:B64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'B6:B64');
                                        % LOG RGDP, UNRATE, INFLATION, FFR
                                        % LOG GOV, LOG RSTPRICES, DLOGRDEBT
                                        
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','B6:B64');
                                       
elseif dtype == 2
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'C6:C64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'B6:B64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','C6:C64');
                                       
elseif dtype == 3
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'D6:D64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'C6:C64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','D6:D64');
                                       
elseif dtype == 4
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'E6:E64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'D6:D64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','E6:E64');
                                       
elseif dtype == 5
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'F6:F64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'E6:E64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','F6:F64');
                                       
elseif dtype == 6
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'G6:G64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'F6:F64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','G6:G64');
                                       
elseif dtype == 7
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'H6:H64');

    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'G6:G64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','H6:H64');
                                       
elseif dtype == 8
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'I6:I64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'H6:H64');
                                       
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','I6:I64');
                                       
else
    
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'J6:J64');
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'I6:I64');
                                      
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','J6:J64');
                                      
end
                                                                           
[data.Var3_Controls,~,~]                = xlsread('DATA_Mertens2015', ...
                                           'CONTROLS','B6:H64');
                                        % LOG RGDP, UNRATE, INFLATION, FFR
                                        % LOG GOV, LOG RSTPRICES, DLOGRDEBT
                                                                             
cd ..

clear ctype dtype ok1 ok2 prompt1 prompt2 strs1 strs2

%% 3) Least-squares, reduced-form estimation

disp('-')

disp('3) The third section estimates the reduced-form VAR parameters');

disp('(output saved in "RForm" structure)')

addpath(strcat(main_d,'/functions/RForm'));

SVARinp.ydata    = ...
    [-log(1-data.Var1_AMTR),data.Var2_LogIncome,data.Var3_Controls];

SVARinp.Z        = data.Var4_ExtIV;

SVARinp.n        = size(SVARinp.ydata,2);

RForm.p          = p; %RForm.p is the number of lags in the model


%a) Estimation of (AL, Sigma) and the reduced-form innovations

[RForm.mu, ...
 RForm.AL, ...
 RForm.Sigma,...
 RForm.eta,...
 RForm.X,...
 RForm.Y]        = RForm_VAR(SVARinp.ydata,p);

%b) Estimation of Gammahat (n times 1)

RForm.Gamma      = RForm.eta*SVARinp.Z(p+1:end,1)/(size(RForm.eta,2)); 
%(We need to take the instrument starting at period (p+1), because
%we there are no reduced-form errors for the first p entries of Y.)

%c) Add initial conditions and the external IV to the RForm structure
    
RForm.Y0         = SVARinp.ydata(1:p,:);
    
RForm.externalIV = SVARinp.Z(p+1:end,1);
    
RForm.n          = SVARinp.n;
    
%d) Definitions for next section

    n            = RForm.n;
    
    T            = (size(RForm.eta,2));
    
    k            = size(RForm.Gamma,2);
    
    d            = ((n^2)*p)+(n);     %This is the size of (vec(A)',Gamma')'
    
    dall         = d+ (n*(n+1))/2;    %This is the size of (vec(A)',vech(Sigma), Gamma')'
    
display(strcat('(total number of parameters estimated:',num2str(d),'; sample size:',num2str(T),')'));

%% 4) Estimation of the asymptotic variance of A,Gamma

disp('-')

disp('4) The fourth section estimates the asymptotic covariance matrix of the reduced-form VAR parameters');

disp('(output saved in "RForm" structure)')

%a) Covariance matrix for vec(A,Gammahat). Used
%to conduct frequentist inference about the IRFs. 

[RForm.WHatall,RForm.WHat,RForm.V] = ...
    CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,nw);                

%NOTES:
%The matrix RForm.WHatall is the covariance matrix of 
% vec(Ahat)',vech(Sigmahat)',Gamma')'
 
%The matrix RForm.WHat is the covariance matrix of only
% vec(Ahat)',Gamma')' 
 
% The latter is all we need to conduct inference about the IRFs,
% but the former is needed to conduct inference about FEVDs. 

%% 5) Compute standard and weak-IV robust confidence set suggested in MSW

disp('-')

disp('5) The fifth section reports standard and weak-IV robust confidence sets ');

disp('(output saved in the "Inference.MSW" structure)')

%Set-up the inputs for the MSW function

nvar   =  1;  %Variable used for normalization

x      = -1;  %Scale of the shock

hori   =  6;  %Number of horizons for the IRFs

%Apply the MSW function

tic;

addpath(strcat(main_d,'/functions'));

[InferenceMSW,Plugin,Chol] = MSWfunction(confidence,nvar,x,hori,RForm,1);

disp('The MSW routine takes only:')

toc;

%Report the estimated shock:
    
epsilonhat=Plugin.epsilonhat;

epsilonhatstd=Plugin.epsilonhatstd;

%% 6) Plot Results

addpath(strcat(main_d,'/functions/figuresfun'));

figure(1)

plots.name(1,:)  = {'Log(1/1-AMTR)'};

plots.name(2,:)  = {'Log Income'};

plots.name(3,:)  = {'Log Real GDP'};

plots.name(4,:)  = {'Unemployment Rate'};

plots.order      = [1,3,2,4];

caux             = norminv(1-((1-confidence)/2),0,1);

for iplot = 1:4
    
    subplot(2,2,plots.order(1,iplot));
    
    plot(Plugin.IRF(iplot,:),'b'); hold on
    
    [~,~] = jbfill(1:1:hori+1,InferenceMSW.MSWubound(iplot,:),...
        InferenceMSW.MSWlbound(iplot,:),[204/255 204/255 204/255],...
        [204/255 204/255 204/255],0,0.5); hold on
    
    dmub  =  Plugin.IRF(iplot,:) + (caux*Plugin.IRFstderror(iplot,:));
    
    lmub  =  Plugin.IRF(iplot,:) - (caux*Plugin.IRFstderror(iplot,:));
    
    h1 = plot(dmub,'--b'); hold on
    
    h2 = plot(lmub,'--b'); hold on
    
    clear dmub lmub
    
    h3 = plot([1 6],[0 0],'black'); hold off
    
    xlabel('Year')
    
    title(plots.name(iplot,:));
    
    if iplot == 1
        
        legend('SVAR-IV Estimator',strcat('MSW C.I (',num2str(100*confidence),'%)'),...
            'D-Method C.I.')
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','southeast')
     
    end
      
    if iplot == 1
        
        axis([1 6 -1.4 .6]);
        
    elseif iplot == 2
            
        axis([1 6 -.5 2.5]);
            
    elseif iplot == 3
        
        axis([1 6 -.4 1.6]);
                
    elseif iplot == 4
                    
        axis([1 6 -.7 .3]);
    end
    
end

%% 7) Save the output and plots in ./Output/Mat and ./Output/Figs

disp('-')

disp('7) The final section saves the .mat files and figures in the Output folder')

%Check if the Output File exists, and if not create one.

if exist('Output','dir')==0
    
    mkdir('Output')
        
end

if exist('Output/Mat','dir')==0
    
    mkdir('Output/Mat')
        
end

if exist('Output/Figs','dir')==0
    
    mkdir('Output/Figs')
        
end

cd(strcat(main_d,'/Output/Mat'));

output_label = strcat(date,'_p=',num2str(p));

save(strcat('karel_68_all_1950_2008_',output_label,'.mat'),...
     'InferenceMSW','Plugin','RForm','SVARinp');

figure(1)
 
cd(strcat(main_d,'/Output/Figs'));

print(gcf,'-depsc2',strcat('karel_68_all_1950_2008_',output_label,'.eps'));

cd(main_d);

clear plots output_label main_d k e_LP
