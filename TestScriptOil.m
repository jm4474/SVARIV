%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This version: August 28th, 2018
% Comment: We have tested this function on an iMac 
%         @3.4 GHz Intel Core i5 (16 GB 2400 MHz DDR4)
%         Running Matlab R2018a.
%         This script runs in about 10 seconds. 

clear;

tic

direct = pwd;

%% 1) Set number of VAR Lags, Newey West lags and confidence level.

fprintf('This script is an example on how to use the function SVARIV_General to report confidence intervals for  IRFs \n');

fprintf('estimated using the SVAR-IV approach described in MSW(18).\n');

disp('-');

disp('Section 1 describes sets values for VAR Lags, Newey West lags (if desired), confidence level, VAR variable names, time unit,');

disp('variable used for normalization, scale of the shock and the number of horizons for the IRFs. It also allows the user');

disp('to select specific variables of interest for the IRF plots and cumulative IRF plots');

application = 'Oil';  % Name of this empirical application. This name will be used for creating and accessing folders

p           = 24;     %Number of lags in the VAR model
 
confidence  = .95;    %Confidence Level for the standard and weak-IV robust confidence set

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
                 
IRFselect   = [1];
% By default, the program generates a single figure with the IRFs for ALL variables
% in the VAR. However, IRFselect allows the user to generate an indepedent
% figure displaying only some specific variables of interest. 
% The program also saves the IRFs in IRFselect as separate .eps files 

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the above vector will select the variables 
% "Index of real economic activity" and the "Real Price of Oil".

cumselect = [2]; 
%cumselect allows the user to generate cumulative IRF plots displaying only 
%some specific variables of interest. 

% Make sure to match the indices above to the variables in your
% dataset. E.g. the above vector will select the variable "Percent Change 
% in Global Crude Oil Production."

%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

disp('-')

disp('Section 2 in this script loads the data necessary for the SVAR-IV analysis and specifies the name of the dataset.')


cd(strcat(pwd,'/Data/',application));

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

cd(direct)
 
%% 3) Test

disp('-')

disp('Section 3 in this script calls the SVARIV_General function')

savdir = strcat(direct,'/Output/',application);  %selected directory where the output files will be saved

addpath(strcat(direct,'/functions/MasterFunction')); 

addpath(strcat(direct,'/functions/Inference'));

[Plugin, InferenceMSW, Chol, RForm, figureorder, grid] = SVARIV_General(p, confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name);
 
% A more in depth description of the function can be found within the
% function file. For clarity purposes, we briefly describe the function
% below:
%
% -Inputs:
%       p:            Number of lags in the VAR model                                                    (1 times 1)                                          
%       confidence:   Value for the standard and weak-IV robust confidence set                           (1 times 1) 
%       ydata:        Endogenous variables from the VAR model                                            (T times n) 
%       z:            External instrumental variable                                                     (T times 1)
%       NWlags:       Newey-West lags                                                                    (1 times 1)
%       norm:         Variable used for normalization                                                    (1 times 1)
%       scale:        Scale of the shock                                                                 (1 times 1)
%       horizons:     Number of horizons for the Impulse Response Functions(IRFs)                        (1 times 1)
%       savdir:       Directory where the figures generated will be saved                                (String)
%       columnnames:  Vector with the names for the endogenous variables, in the same order as ydata     (1 times n)
%       IRFselect:    Indices for the variables that the user wants separate IRF plots for               (1 times q)
%       time:         Time unit for the dataset (e.g. year, month, etc.)                                 (String)
%       dataset_name: The name of the dataset used for generating the figures (used in the output label) (String)
%
% -Output:
%       Plugin: Structure containing standard plug-in inference
%       InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%       Chol: Cholesky IRFs

%% 4) "Standard" bootstrap-type inference based on samples from the asy dist.

disp('Section 4 in this script calls the Gasydistboots function to provide inference for SVAR-IV based on samples from the asy. dist.')

seed            = load(strcat(direct,'/seed/seedMay12.mat')); 
    
seed            = seed.seed;

n = RForm.n; 

T  = size(ydata, 1);            % Number of observations/time periods.

NB = 1000;                      % Number of bootstrap replications 

SVARinp.Z = z;

SVARinp.ydata = ydata;

SVARinp.n = n;


addpath('functions/Inference');

[~,InferenceMSW.bootsIRFs] = ...
                  Gasydistboots(seed, NB, n, p, norm, scale, horizons, confidence, T,...
                  @IRFSVARIV, SVARinp, NWlags, RForm.AL, RForm.Sigma, RForm.Gamma, RForm.V, RForm.WHatall);

%% 5) Bootstrap Plots

disp('Section 5 in this script calls the Bootstrap_Plots function.')

addpath(strcat(direct,'/functions/figuresfun'));

[caux,InferenceMSW] = Bootstrap_Plots(n, p, horizons, confidence, RForm, SVARinp, figureorder, Plugin, InferenceMSW, time, columnnames, savdir, direct, dataset_name, IRFselect, cumselect);


%% 6) AR confidence set using bootstrap implementation

disp('Section 6 in this script calls the GasydistbootsAR function to do the bootstrap implementation of the Anderson-Rubin confidence set')

cd(strcat(direct,'/functions/Inference'));

[reject, bootsIRFs] = GasydistbootsAR(ydata, T, seed, RForm.n, NB, p, norm, scale, horizons, confidence, SVARinp, NWlags, RForm.AL, RForm.Sigma, RForm.Gamma, RForm.V, RForm.WHatall, grid);

toc; 
