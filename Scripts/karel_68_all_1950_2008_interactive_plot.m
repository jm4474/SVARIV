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

%% 3) Least-squares, reduced-form estimation

disp('-')

disp('3) The third section estimates the reduced-form VAR parameters');

disp('(output saved in "RFform" structure)')

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

%% 6) Compute the MSW confidence set
%------------------------------------------
%(output saved in the "Inference.MSW" structure)
%------------------------------------------

tic;

disp('d) Some comments regarding the MSW inference procedure:')

%Set-up the inputs for the MSW function

nvar   = 1;

x      = -1;

hori   = 6;

%Apply the MSW function

addpath(strcat(main_d,'/functions'));

[InferenceMSW,Plugin,Chol] = MSWfunction(confidence,nvar,x,hori,RForm,1);

%Report the estimated shock:
    
%epsilonhat=Plugin.epsilonhat;

%epsilonhatstd=Plugin.epsilonhatstd;

disp('The MSW routine takes only:')

toc;

%% Extra

%% Draws from the reduced-form parameters to conduct "Standard" Bootstrap inference
%------------------------------------------
%(output saved in the "Inference" structure)
%------------------------------------------

    %a) Generate Samples from vec(A,Gammahat)
    
    seed = load(strcat(main_d,'/seed/seedMay12.mat')); %The seed file is in the Seed folder
    
    seed = seed.seed;
    
    rng(seed); clear seed
    
    %Inference.gvar=[randn(dall,1000),zeros(dall,1)];
    
    %Adjust the covariance matrix to make it symmetric
    
    RForm.Whatalladj  = (RForm.WHatall+RForm.WHatall')/2;
    
    [aux1,aux2]       = eig(RForm.Whatalladj);
    
    RForm.Whatallpos  = aux1*max(aux2,0)*aux1';
    
    Inference.gvar    = [mvnrnd(zeros(1000,dall),(RForm.Whatallpos)/T)',zeros(dall,1)];
    
    Inference.Draws   = bsxfun(@plus,Inference.gvar,...
                        [RForm.AL(:);RForm.V*RForm.Sigma(:);RForm.Gamma(:)]);

    %The vector "Inference.Draws" represents a vector of 1000
    %from a multivariate normal vector (of dimension dall) centered
    %at [vec(A)', vech(Sigma)',Gamma']' with covariance matrix 
    %(RForm.WHatall/T). Thus, it represents a draw from the asy. dist
    %of the reduced-form parameters.    
    
    
%% 8) Map from (A,Gamma) to IRF and from (A,Sigma,Gamma to IRFs)
 
%   Definitions for this section

     I = size(Inference.Draws,2);
     
     disp('f) The number of draws from the asymptotic distribution of vec(A,Gamma)')
     
     display(strcat('used to implement bootstrap is:',num2str(I)));

    %a) Evaluate the IRFs at each draw of A,Gamma
    e=eye(n); %the columns of this identity matrix will be used to evaluate IRFs and FEVDs 
    
    tic;
    
    disp('The numerical procedure to implement standard inference is running...')  
      
   for ip = 1:I
      %i) Generate the draws for AL and check that they fall in the 
      %   stationarity region
      
      AL        = reshape(Inference.Draws(1:(n^2)*p,ip),[n,n*p]);
      
      %ii) Generate the draws from Sigma and check they are positive definite
      
      vechSigma = Inference.Draws((n^2)*p+1:(n^2)*p+(n*(n+1)/2),ip);
      
      Sigma     = tril(ones(n),0); Sigma(Sigma==1) = vechSigma';
      
      Sigma     = Sigma + tril(Sigma,-1)';
      
      %This is a simple way to create a matrix Sigma from the matrix vech(Sigma)
      
      if min(eig(Sigma))>0
          
          Inference.pdSigma(:,ip) = 1;
          
      else
          
          Inference.pdSigma(:,ip) = 0;
          
      end
      
      %iii) Draws from Gamma
      
      Gamma     = reshape(Inference.Draws(((n^2)*p)+(n*(n+1)/2)+1:end,ip),[n,k]);
     
      %iV) Reduced-form MA coefficients
      
      Cauxsim   = [eye(n),MARep(AL,p,hori)]; 
      
      Csim      = reshape(Cauxsim,[n,n,hori+1]);
      
      
      %V)  Obtain the plug-in estimator for each draw
      
      B1        = x*Gamma./Gamma(nvar,1);    
      
      %v)   Generate the normalized IRFs for each draw
      %1 is for standard IRF. 2 is for cumulative
      
      Inference.IRFZ(:,:,ip,1) = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);
      
      clear AL vechSigma Gamma Csim Cauxsim  B1 B1aux
      
   end
   
    disp('The numerical procedure is over.')
    
    toc;

%% 9) Implement "Standard" Bootstrap Inference

   aux                  = reshape(Inference.pdSigma,[1,1,I]);
   
   Inference.bootsIRFZ  = quantile(Inference.IRFZ(:,:,aux==1,:),...
                          [((1-confidence)/2),1-((1-confidence)/2)],3);
     
%% 10)

addpath(strcat(main_d,'/functions/figuresfun'));


%% Plots

figure(1)

plots.name(1,:)  = {'Log(1/1-AMTR)'};

plots.name(2,:)  = {'Log Income'};

plots.name(3,:)  = {'Log Real GDP'};

plots.name(4,:)  = {'Unemployment Rate'};

plots.order       = [1,3,2,4];

for iplot = 1:4
    
    subplot(2,2,plots.order(1,iplot));
    
    plot(Plugin.IRF(iplot,:),'b'); hold on
    
    [~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(iplot,:,2,1),...
        Inference.bootsIRFZ(iplot,:,1,1),[204/255 204/255 204/255],...
        [204/255 204/255 204/255],0,0.5); hold on
    
    h1 = plot(InferenceMSW.MSWubound(iplot,:),'--b'); hold on
    
    h2 = plot(InferenceMSW.MSWlbound(iplot,:),'--b'); hold on
    
    h3 = plot([1 6],[0 0],'black'); hold off
    
    xlabel('Year')
    
    title(plots.name(iplot,:));
    
    if iplot == 1
        
        legend('Point Estimator','Asy Dist. Boots','MSW-WeakIV')
        
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        legend boxoff
        
        legend('location','northwest')
     
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

%%

if exist('Output','dir')==0
    
    mkdir('Output')
        
end

cd(strcat(main_d,'/Output'));

print(gcf,'-depsc2','karel_68_all_1950_2008_interactive.eps');


%% 


save('karel_68_all_1950_2008_interactive.mat');






