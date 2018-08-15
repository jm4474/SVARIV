%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This is a test version to see whether the SVARIV_Luigi function is working: Auguest 13th, 2018
% Comment: We have tested this script on a Macbook Pro 
 
direct = pwd;

%% 1) Set number of VAR Lags, Newey West lags and confidence level.
 
p           = 2; 
 
confidence  = .68;

% Define the variables in the SVAR
columnnames = [{'Log(1/1-AMTR)'}, ...
               {'Log Income'}, ...
               {'Log Real GDP'},...
               {'Unemployment Rate'},...
               {'a'}, {'b'}, {'c'}, {'d'}, {'e'}];

time        = 'Year';

NWlags      = 8;

norm        = 1;

scale       = -1;

horizons    = 6;
 
%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

cd(strcat(pwd,'/Data'));
 
    [data.years, ~, ~]                  = xlsread('DATA_Mertens2015',... 
                                          'AMTR (Figure 1)', 'A6:A64');
                                        % time variable. Here as a 59x1
                                        % vector
 
    [data.Var1_AMTR, ~, ~]              = xlsread('DATA_Mertens2015',...
                                          'AMTR (Figure 1)', 'C6:C64');
                                        % 59x1 vector.
                                      
    [data.Var2_LogIncome,~,~]           = xlsread('DATA_Mertens2015',...
                                          'LOG AVG INCOME', 'B6:B64');
                                        % 59x1 vector
                                      
    [data.Var3_Controls,~,~]            = xlsread('DATA_Mertens2015', ...
                                           'CONTROLS','B6:H64');
                                        % 59x7 matrix with the control
                                        % variables. Here they are:
                                        % LOG RGDP, UNRATE, INFLATION, FFR
                                        % LOG GOV, LOG RSTPRICES, DLOGRDEBT
                                        
    [data.Var4_ExtIV,~,~]               = xlsread('DATA_Mertens2015',...
                                          'PROXIES (Table 3)','C6:C64');
                                        %External IV. 59x1 vector.
                                      
cd ..
 
%% 3 Test
 
ydata    = ...
    [-log(1-data.Var1_AMTR),data.Var2_LogIncome,data.Var3_Controls];
 
z        = data.Var4_ExtIV;
 

savdir = strcat(direct,'/Output');
 
addpath(strcat(direct,'/functions'));

[Plugin, InferenceMSW, Chol] = SVARIV_General(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, time);
 
%To Do: Provide a description of each of the inputs, so that the user knows
%       what are we talking about. 
