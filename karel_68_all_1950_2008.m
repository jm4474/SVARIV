%% 1) READ ME

% Companion Script: Marginal Tax Rates and Income: New Time Series evidence. 

% We suggest you to run this script section by section. 
beep off
clear; clc;

% Change the main directory if needed

main_d = pwd;

cd(main_d);

% OPTIONS:
% The following options are avaiable for this script:    
% -confidence level:
confidence= .68;
% -Name of the external instrument:
instrumentname='Karel';
%This program was tested using Matlab R2015a and a Macbook Pro
%2.4 GHz Intel Core i7 with OS X El Capitan.

%Before running this program make sure that this .m file is in the same
%directory as the folders:

%AuxFunctions
%Data
%Seed

%This m.file will not run properly without the input in the folders above.

%% 2) Loading the data 
%--------------------------------------
%(Ouput saved in the "dataU" structure)
%(NOTE: we decided to use "structures" just to keep the 
%workspace organized)
%-------------------------------------- 
cd(strcat(main_d,'/Data'));

%a) All units
[data.years,~,~]=xlsread('DATA_Mertens2015','AMTR (Figure 1)','A6:A64');
[data.Var1_AMTR,~,~]=xlsread('DATA_Mertens2015','AMTR (Figure 1)','C6:C64');
[data.Var2_LogIncome,~,~]=xlsread('DATA_Mertens2015','LOG AVG INCOME','B6:B64');
[data.Var3_Controls,data.Var3_Labels,~]=xlsread('DATA_Mertens2015','CONTROLS','B6:H64');
[data.Var4_ExtIV,~,~]=xlsread('DATA_Mertens2015','PROXIES (Table 3)','C6:C64');

cd ..

%% 3) Least-squares, reduced-form estimation

addpath(strcat(main_d,'/3Auxfunctions/RForm'));
%SVARinp.ydata =[log(1-data.Var1_AMTR),data.Var2_LogIncome,data.Var3_Controls];
SVARinp.ydata =[-log(1-data.Var1_AMTR),data.Var2_LogIncome,data.Var3_Controls];
SVARinp.Z = data.Var4_ExtIV;

SVARinp.n = size(SVARinp.ydata,2);
p=2;

%% 4) Reduced-form VAR estimation 
%---------------------------------------
%(Output saved in the "RForm" Structure)    
%---------------------------------------

RForm.p=p; 


%a) Estimation of (AL, Sigma) and the reduced-form innovations
    %RForm.p is the number of lags in the model
    %NOTE: You can replace the function "RForm_VAR.m" by your own Matlab
    %function to estimate reduced-form parameters. 
      
    [RForm.mu,RForm.AL,RForm.Sigma,RForm.eta,RForm.X,RForm.Y] = RForm_VAR(SVARinp.ydata,p);

%b) Estimation of Gammahat
    RForm.Gamma= RForm.eta*SVARinp.Z(p+1:end,1)/(size(RForm.eta,2)); %r times k
    %We need to take the instrument starting at period (p+1), because
    %we only report the reduced-form errors for the first p entries of Y.

%c) Add initial conditions and the external IV to the RForm structure
    
    RForm.Y0=SVARinp.ydata(1:p,:);
    RForm.externalIV=SVARinp.Z(p+1:end,1);
    RForm.n=SVARinp.n;
    
%d) Definitions for next section
    n=RForm.n; T=(size(RForm.eta,2)); k=size(RForm.Gamma,2);
    d=((n^2)*p)+(n);     %This is the size of (vec(A)',Gamma')'
    dall=d+ (n*(n+1))/2; %This is the size of (vec(A)',vech(Sigma), Gamma')'
    display(strcat('c) Total number of parameters estimated:',num2str(d)));
    display(strcat('(and the sample size is T=',num2str(T),')'));  

%% 5) Estimation of the asymptotic variance of A,Gamma


%a) Covariance matrix for vec(A,Sigma,Gammahat). This matrix will be used
%to conduct frequentist inference about the IRFs. 

    [RForm.WHatall,RForm.WHat,RForm.V] = CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,8);                
  
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
display('d) Some comments regarding the MSW inference procedure:')

%Set-up the inputs for the MSW function
nvar=1;
x=-1;
hori=6;

%Apply the MSW function

addpath(main_d);
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
    seed=load(strcat(main_d,'/4Seed/seedMay12.mat')); %The seed file is in the Seed folder
    seed=seed.seed;
    rng(seed); clear seed
    %Inference.gvar=[randn(dall,1000),zeros(dall,1)]; 
    %Adjust the covariance matrix to make it symmetric
    RForm.Whatalladj=(RForm.WHatall+RForm.WHatall')/2;
    [aux1,aux2]=eig(RForm.Whatalladj);
    RForm.Whatallpos=aux1*max(aux2,0)*aux1';
    Inference.gvar=[mvnrnd(zeros(1000,dall),(RForm.Whatallpos)/T)',zeros(dall,1)];  
    Inference.Draws=...
    bsxfun(@plus,Inference.gvar,...
    [RForm.AL(:);RForm.V*RForm.Sigma(:);RForm.Gamma(:)]);

    %The vector "Inference.Draws" represents a vector of 1000
    %from a multivariate normal vector (of dimension dall) centered
    %at [vec(A)', vech(Sigma)',Gamma']' with covariance matrix 
    %(RForm.WHatall/T). Thus, it represents a draw from the asy. dist
    %of the reduced-form parameters.    
    
    
%% 8) Map from (A,Gamma) to IRF and from (A,Sigma,Gamma to IRFs)
 
%   Definitions for this section
     I=size(Inference.Draws,2);
     disp('f) The number of draws from the asymptotic distribution of vec(A,Gamma)');
     display(strcat('used to implement bootstrap is:',num2str(I)));

    %a) Evaluate the IRFs at each draw of A,Gamma
    e=eye(n); %the columns of this identity matrix will be used to evaluate IRFs and FEVDs     
    tic
    disp('The numerical procedure to implement standard inference is running...')  
      
   for ip=1:I
      %i) Generate the draws for AL and check that they fall in the 
      %   stationarity region
      AL=reshape(Inference.Draws(1:(n^2)*p,ip),[n,n*p]);      
      %ii) Generate the draws from Sigma and check they are positive
      %    definite
      vechSigma=Inference.Draws((n^2)*p+1:(n^2)*p+(n*(n+1)/2),ip);
      Sigma=tril(ones(n),0); Sigma(Sigma==1)=vechSigma';
      Sigma=Sigma + tril(Sigma,-1)'; 
      %This is a simple way to create a matrix Sigma from the matrix 
      %vechSigma
      if min(eig(Sigma))>0
      Inference.pdSigma(:,ip)=1;
      else
      Inference.pdSigma(:,ip)=0;
      end
      
      %iii) Draws from Gamma
      Gamma=reshape(Inference.Draws(((n^2)*p)+(n*(n+1)/2)+1:end,ip),[n,k]);
     
      %iV) Reduced-form MA coefficients
      Cauxsim=[eye(n),MARep(AL,p,hori)]; 
      Csim= reshape(Cauxsim,[n,n,hori+1]);
      
      
      %V)  Obtain the plug-in estimator for each draw
      
      B1=x*Gamma./Gamma(nvar,1);    
      
      %v)   Generate the normalized IRFs for each draw
      %1 is for standard IRF. 2 is for cumulative 
      Inference.IRFZ(:,:,ip,1)=reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);
      
      
      clear AL vechSigma Gamma Csim Cauxsim  B1 B1aux
   end
    disp('The numerical procedure is over.')
    toc   

%% 9) Implement "Standard" Bootstrap Inference

   aux=reshape(Inference.pdSigma,[1,1,I]);
   %Inference.bootsIRFZ=quantile(Inference.IRFZ(:,:,aux==1,:),[.16,.84],3);
   Inference.bootsIRFZ=quantile(Inference.IRFZ(:,:,aux==1,:),[((1-confidence)/2),1-((1-confidence)/2)],3);
     
%% 10) 
addpath(strcat(main_d,'/3Auxfunctions/figuresfun'));


%% Plots

figure(1)
subplot(2,2,1)
plot(Plugin.IRF(1,:),'b'); hold on
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(1,:,2,1),Inference.bootsIRFZ(1,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
plot(InferenceMSW.MSWubound(1,:),'--b'); hold on
h1=plot(InferenceMSW.MSWlbound(1,:),'--b'); hold on
h2=plot([1 6],[0 0],'black'); hold off
xlabel('Year')
title('1/(1-AMTR)')
legend('Point Estimator','Asy Dist. Boots','MSW-WeakIV')
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend boxoff
legend('location','northwest')
axis([1 6 -1.4 .6])

subplot(2,2,2)
plot(Plugin.IRF(3,:),'b'); hold on
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(3,:,2,1),Inference.bootsIRFZ(3,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
plot(InferenceMSW.MSWubound(3,:),'--b'); hold on
plot(InferenceMSW.MSWlbound(3,:),'--b'); hold on
plot([1 6],[0 0],'black');
hold off
xlabel('Year')
title('Log Real GDP')
axis([1 6 -.4 1.6])

subplot(2,2,3)
plot(Plugin.IRF(2,:),'b'); hold on 
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(2,:,2,1),Inference.bootsIRFZ(2,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
plot(InferenceMSW.MSWubound(2,:),'--b'); hold on
plot(InferenceMSW.MSWlbound(2,:),'--b'); hold on
plot([1 6],[0 0],'black'); hold off
xlabel('Year')
title('Log Income')
axis([1 6 -.5 2.5])

subplot(2,2,4)
plot(Plugin.IRF(4,:),'b'); hold on
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(4,:,2,1),Inference.bootsIRFZ(4,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
plot(InferenceMSW.MSWubound(4,:),'--b'); hold on
plot(InferenceMSW.MSWlbound(4,:),'--b'); hold on
plot([1 6],[0 0],'black'); hold off
xlabel('Year')
title('Unemployment Rate')
axis([1 6 -.7 .3])


%%
cd(strcat(main_d,'/Output'));

print(gcf,'-depsc2','karel_68_all_1950_2008.eps');


%% 


save('karel_68_all_1950_2008.mat');






