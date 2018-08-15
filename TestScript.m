%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This is a test version to see whether the SVARIV_Luigi function is working: Auguest 13th, 2018
% Comment: We have tested this script on a Macbook Pro 
 
%% 1) Set number of VAR Lags, Newey West lags and confidence level.
 
p = 2; 
 
confidence = .68;

columnnames = ["Log(1/1-AMTR)", "Log Income", "Log Real GDP", "Unemployment Rate", "a", "b", "c", "d", "e"];

time = 'Year';

NWlags = 8;

norm = 1;

scale = -1;

horizons = 6;
 
%% 2) Load data (saved in structure "data")

cd('/Users/luigicaloi/Documents/Pepe/SVARIV/Data');
 
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
 
disp('3) The third section estimates the reduced-form VAR parameters');
 
ydata    = ...
    [-log(1-data.Var1_AMTR),data.Var2_LogIncome,data.Var3_Controls];
 
z        = data.Var4_ExtIV;
 

savdir = '/Users/luigicaloi/Documents/Pepe/SVARIV/Output';
 
cd '/Users/luigicaloi/Documents';

addpath('/Users/luigicaloi/Documents/Pepe/SVARIV/functions');

[Plugin, InferenceMSW, Chol] = SVARIV_Luigi(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, time);
 

