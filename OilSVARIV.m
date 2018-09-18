%% This script file implements standard and weak-IV robust SVAR-IV inference.
% This version: August 28th, 2018
% Comment: We have tested this function on an iMac 
%         @3.4 GHz Intel Core i5 (16 GB 2400 MHz DDR4)
%         Running Matlab R2018a.
%         This script runs in about 12 seconds. 

clear;

tic

direct = pwd;

disp('(We would like to thank Luigi Caloi, Qifan Han, Hamza Husain and Jianing Zhai for excellent research assistance)')

%% 1) Set number of VAR Lags, Newey West lags and confidence level.

fprintf('\nThis script uses the function SVARIV to report confidence intervals for  IRFs \n');

fprintf('as described in Montiel Olea, Stock, and Watson(18).\n');

disp('--');

fprintf('Section 1 sets values for: \n')

fprintf('-VAR Lags \n')
fprintf('-Newey West lags (if desired) \n')
fprintf('-Confidence level \n')
fprintf('-VAR variable names \n')
fprintf('-Time unit \n')
fprintf('-Variable used for normalization \n')
fprintf('-Scale of the shock \n')
fprintf('-Number of horizons for the IRFs. \n')

disp('This section also allows the user to select specific variables of interest for the IRF plots');

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
                 
IRFselect   = [1,2];
% By default, the program generates a single figure with the IRFs for ALL variables
% in the VAR. However, IRFselect allows the user to generate an indepedent
% figure displaying only some specific variables of interest. 
% For example, IRFselect   = [1,2] generates individual IRF figures for
% variables 1 and 2 and saves them as separate .eps files 

% Make sure to match the indices above and to the variables in your
% dataset. E.g. the IRFselect   = [1,2] will select the variables 
% "Index of real economic activity" and the "Real Price of Oil".

cumselect = [2,3]; 
%cumselect allows the user to generate cumulative IRF plots displaying only 
%some specific variables of interest. 

% Make sure to match the indices above to the variables in your
% dataset.

%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

disp('-')

disp('Section 2 loads the data and specifies the name of the dataset.')


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
 
%% 3) Run the SVARIV function

disp('--')

disp('Section 3 in this script calls the SVARIV function')

savdir = strcat(direct,'/Output/',application);  %selected directory where the output files will be saved

addpath(strcat(direct,'/functions/MasterFunction')); 

addpath(strcat(direct,'/functions/Inference'));

[Plugin, InferenceMSW, Chol, RForm, figureorder, ARxlim, ARylim] = SVARIV(p, confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name);
 
% A more in depth description of the function can be found within the
% function file. For clarity purposes, we briefly describe the function
% below:
%
% -Inputs:
%       p:            Number of lags in the VAR model                                                    (1 x 1)                                          
%       confidence:   Value for the standard and weak-IV robust confidence set                           (1 x 1) 
%       ydata:        Endogenous variables from the VAR model                                            (T x n) 
%       z:            External instrumental variable                                                     (T x 1)
%       NWlags:       Newey-West lags                                                                    (1 x 1)
%       norm:         Variable used for normalization                                                    (1 x 1)
%       scale:        Scale of the shock                                                                 (1 x 1)
%       horizons:     Number of horizons for the Impulse Response Functions(IRFs)                        (1 x 1)
%       savdir:       Directory where the figures generated will be saved                                (String)
%       columnnames:  Vector with the names for the endogenous variables, in the same order as ydata     (1 x n)
%       IRFselect:    Indices for the variables that the user wants separate IRF plots for               (1 x q)
%       cumselect:    Indices for the variables that the user wants cumulative IRF plots for             (1 x q)
%       time:         Time unit for the dataset (e.g. year, month, etc.)                                 (String)
%       dataset_name: The name of the dataset used for generating the figures (used in the output label) (String)
%
% -Output:
%       Plugin: Structure containing standard plug-in inference
%       InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%       Chol: Cholesky IRFs

%% 4) "Standard" bootstrap-type inference based on samples from the asy dist.

disp('Section 4 calls the Gasydistboots function to provide inference for SVAR-IV based on samples from the asy. dist.')

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
                  @IRFSVARIV, RForm.AL, RForm.Sigma, RForm.Gamma, RForm.V, RForm.WHatall,SVARinp, NWlags);

%% 5) Bootstrap Plots

disp('Section 5 calls the Bootstrap_Plots function.')

addpath(strcat(direct,'/functions/figuresfun'));

[caux,InferenceMSW, figureorder] = Bootstrap_Plots(n, p, horizons, confidence, RForm, SVARinp, figureorder, Plugin, InferenceMSW, time, columnnames, savdir, direct, dataset_name, IRFselect, cumselect);


%% 6) AR confidence set using bootstrap implementation

disp('Section 6 in this script calls the GasydistbootsAR function to do the bootstrap implementation of the Anderson-Rubin confidence set')

cd(strcat(direct,'/functions/Inference'));

multiplier = 1.5;     %Scalar that GasydistbootsAR will use to create the "grid" of null hypothesis for IRFs
                     %IRFhat +- multiplier*ARylim 

grid_size = 50;      % Number of mull hypotheses (lambdas) in the grid, for each variable and for each horizon.

[reject, bootsIRFs, gridpointupperMSW, gridpointlowerMSW, null_grid] = GasydistbootsAR(ydata, T, seed, RForm.n, NB, p, norm, scale, horizons, confidence, SVARinp, NWlags, RForm.AL, RForm.Sigma, RForm.Gamma, RForm.V, RForm.WHatall, Plugin, multiplier, grid_size, ARylim);

%% 7) AR Bootstrap plots

disp('Section 7 calls the BootstrapAR_Plots function.')

addpath(strcat(direct,'/functions/figuresfun'));

BootstrapAR_Plots(n, p, horizons, confidence, RForm, SVARinp, figureorder, Plugin, InferenceMSW, time, columnnames, savdir, direct, dataset_name, IRFselect, cumselect, reject, null_grid, norm);

toc; 
