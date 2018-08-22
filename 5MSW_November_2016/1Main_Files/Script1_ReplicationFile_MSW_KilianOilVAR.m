%Replication File. Last update: November 21st, 2016. 
%% 1) READ ME

% This script estimates the three-dimensional oil SVAR proposed by
% Kilian (2009; "Not oil shocks are alike ...") using Kilian's (2008)
% short-fall of OPEC oil production as an external IV.

% To estimate the SVAR this script applies the procedures described in 
%Montiel-Olea, Stock and Watson (2016). 

% We suggest you to run this script section by section. 
beep off
clear; clc;

%The VAR under consideration includes 3 variables: 
%   -Growth Rate of World Oil Production (Source: Kilian 2009)
%   -Global Real Activity Index (Source: Kilian 2009)
%   -Real Price of Oil (Source: Kilian 2009)

% The data can be found here:
% https://www.aeaweb.org/articles?id=10.1257/aer.99.3.1053

% OPTIONS:
% The following options are avaiable for this script:    
% -confidence level:
confidence= .95;
% -Name of the external instrument:
instrumentname='K';
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

main_directory=pwd;
cd 2Data/Kilian_AER_data

data.Kilian=load('data.txt'); 
%The frequency of this data is 1973:2 - 2007:12
%The file data.txt was obtained directly from the AER website

cd(main_directory)
data.ExternalIV=load('2Data/ExternalIV.txt');
%The frequency of this data is 1973:2 - 2004:09
%The .xls file was created by Jim Stock and Mark Watson and it 
%contains a monthly measure of Kilian's [2008] instrument for the
%oil supply shock. 

%a) Name the 3 variables in the VAR
    dataU.varnames = cell(3,1);
    dataU.varnames{1,1} = 'Percent Change in Global Crude Oil Production';
    dataU.varnames{2,1} = 'Index of real economic activity';
    dataU.varnames{3,1} = 'Real Price of Oil';
    
    
%Here we decided to display some code to remind the user what the 
%program is doing.

    display('This program estimates a SVAR using monthly data on:')
    display('1) Percent Change in Global Crude Oil Production');
    display('2) Index of real economic activity'); 
    display('3) Real Price of Oil');  
        
    
%b) Set the Dates: February 1973 to September 2004

    [Year,Month] = meshgrid(1973:2004, 1:12);
    Dates = datenum([Year(:), Month(:), ones(numel(Year),1)]);
    dataU.dates=datestr(Dates(2:end,1)); 
    clear Dates
    years=Year(:);
    months=Month(:);
    dataU.data(:,1)=years(2:end-3,1);  %years
    dataU.data(:,2)=months(2:end-3,1); %months 
    clear years months Year Month
    dataU.namesvar{1,1}='Year';
    dataU.namesvar{2,1}='Month';  
    
    %initial year and month in the sample
    startyear=1973; startmonth=2;
    
    display(strcat('The sample period is:',num2str(startyear),':',num2str(startmonth),'-2004:09.'));
       
    

%c) Create the dataU structure: dataU stands for data "Used"
    dataU.data(:,3) = data.Kilian(1:end-39,1); 
    dataU.data(:,4) = data.Kilian(1:end-39,2);
    dataU.data(:,5) = data.Kilian(1:end-39,3);
    
    
    
    startdata=find(dataU.data(:,1)==startyear & dataU.data(:,2)==startmonth );
     
%% 3) Parameters for the SVAR exercise 
%--------------------------------------
%(output saved in the "SVARinp" structure)
%--------------------------------------  

addpath('3AuxFunctions/RForm');
% The RForm folder contain different .m files created by 
% Matthias Meier (University of Bonn) and Jos? Luis Montiel Olea
%(Columbia University). These files estimate different objects of 
% interest in reduced-form SVAR analysis. 
 
[bic,aic,hqic]=bicaic(dataU.data(startdata:end,3:5),24);
%The function file bicaic.m estimates the number of lags according to
%bic: Bayes Information Criterion
%aic: Akaike Information Criterion
%hqic:Hannan-Quinn Information Criterion

%The inputs for []=bicaic(series,pmax) are 
%a) series: the time series (stored as "data" in the structure "dataU")
%b) pmax: a bound on the number of lags permitted. We set
p = 24; %(as in the specification suggested by Kilian) 


%We now display some of the information obtained from the function bicaic
display(blanks(1));
display('Comments:');
display(strcat('a) Number of lags suggested by the Akaike Information Criterion:',num2str(aic),'; Number of lags used:',num2str(p)));
clear aic;

%a) Use data from February 1973 to September 2004
     
    SVARinp.ydata =dataU.data(startdata:end,3:5); 

%We create the structure "SVARinp" to collect the inputs for SVAR
%an
   
    
%b) Inclue the instrumental variable in SVARinp

    SVARinp.Z = data.ExternalIV;
    display('b) The model is identified using the instrument in Kilian (08)');
    
        
    [SVARinp.n] = size(SVARinp.ydata,2);
    %SVARinp.n is simply the dimension of the SVAR
    
    clear data;     
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
    
    %Stationarity Check
    if max(abs( eig([RForm.AL; eye(SVARinp.n*(RForm.p-1)),...
        zeros(SVARinp.n*(RForm.p-1),SVARinp.n)]) )) >=1;
        error('The processes are not stationary! Check!');
    else
    end
    
    
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

    [RForm.WHatall,RForm.WHat,RForm.V] = CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta);                
  
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
x=-20;
hori=20;

%Apply the MSW function

cd('1Main_Files')

[InferenceMSW,Plugin,Chol] = MSWfunction(confidence,nvar,x,hori,RForm,1);

cd ..

%Report the estimated shock:
    
epsilonhat=Plugin.epsilonhat;    
epsilonhatstd=Plugin.epsilonhatstd;

display('The MSW routine takes only:')
toc;
         
%% 7) Draws from the reduced-form parameters to conduct "Standard" Bootstrap inference
%------------------------------------------
%(output saved in the "Inference" structure)
%------------------------------------------

    %a) Generate Samples from vec(A,Gammahat)
    seed=load('4Seed/seedMay12.mat'); %The seed file is in the Seed folder
    seed=seed.seed;
    rng(seed); clear seed
    Inference.gvar=[randn(dall,1000),zeros(dall,1)];  
    Inference.Draws=...
    bsxfun(@plus,(RForm.WHatall/T)^(.5)*Inference.gvar,...
    [RForm.AL(:);RForm.V*RForm.Sigma(:);RForm.Gamma(:)]);

    %The vector "Inference.Draws" represents a vector of 1000
    %from a multivariate normal vector (of dimension dall) centered
    %at [vec(A)', vech(Sigma)',Gamma']' with covariance matrix 
    %(RForm.WHatall/T). Thus, it represents a draw from the asy. dist
    %of the reduced-form parameters.    
    
    
%% 8) Map from (A,Gamma) to IRF and from (A,Sigma,Gamma to IRFs)
addpath('3AuxFunctions/StructuralIRF')
 
%   Definitions for this section
     I=size(Inference.Draws,2);
     display('f) The number of draws from the asymptotic distribution of vec(A,Gamma)');
     display(strcat('used to implement bootstrap is:',num2str(I)));

    %a) Evaluate the IRFs at each draw of A,Gamma
    e=eye(n); %the columns of this identity matrix will be used to evaluate IRFs and FEVDs     
    tic
    display('The numerical procedure to implement standard inference is running...')  
      
   for ip=1:I;
      %i) Generate the draws for AL and check that they fall in the 
      %   stationarity region
      AL=reshape(Inference.Draws(1:(n^2)*p,ip),[n,n*p]);
      if max(abs( eig([reshape(AL,[n,n*p]); eye(SVARinp.n*(RForm.p-1)),...
        zeros(SVARinp.n*(RForm.p-1),SVARinp.n)]) )) >=1;
        Inference.stat(:,ip)=0;
      else
        Inference.stat(:,ip)=1;
      end
      
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
      Ccumsim=cumsum(Csim,3);
      
      %V)  Obtain the plug-in estimator for each draw
      
      B1=x*Gamma./Gamma(nvar,1);    
      B1aux=Gamma./((Gamma'*Sigma^(-1)*Gamma)^.5);
      
      %v)   Generate the normalized IRFs for each draw
      %1 is for standard IRF. 2 is for cumulative 
      Inference.IRFZ(:,:,ip,1)=reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);
      Inference.IRFZ(:,:,ip,2)=reshape(sum(bsxfun(@times,Ccumsim,B1'),2),[n,hori+1]);
      
      %vi) Generate the FEVDs for each draw (with an inner loop)
      %(This part of the code is a little messy, but we will try to
      %explain it)
      CauxFEV1=zeros(n,n,hori+1);
      CauxFEV2=zeros(n,n,hori+1);
      for ixvar=1:n  %ixvar is just an index for the SVAR time series
          for ishori=1:hori+1 %Ishori is the index for the horizon
              %Now we implement the formula for FEVDs in our paper
              %(as a function of A, Sigma,Gamma)
              CauxFEV1(:,:,ishori)=Csim(:,:,ishori)'*e(:,ixvar)*e(:,ixvar)'*Csim(:,:,ishori);
              CauxFEV2(:,:,ishori)=e(:,ixvar)'*Csim(:,:,ishori)*Sigma*Csim(:,:,ishori)'*e(:,ixvar);
              auxFEVD=(Sigma^(1/2))*(sum(CauxFEV1(:,:,1:ishori),3)./sum(CauxFEV2(:,:,1:ishori),3))*(Sigma^(1/2));
              Inference.FEVD(ixvar,ishori,ip,1)=max(eig(auxFEVD));
              Inference.FEVD(ixvar,ishori,ip,2)=min(eig(auxFEVD));
              Inference.PointEstFEVD(ixvar,ishori,ip)=B1aux'*(Sigma^(-1/2))*auxFEVD*(Sigma^(-1/2))*B1aux;
              B1cholsim=chol(Sigma)';
              B1cholsim=B1cholsim(:,1);
              Inference.CholFEVD(ixvar,ishori,ip)=B1cholsim'*(Sigma^(-1/2))*auxFEVD*(Sigma^(-1/2))*B1cholsim;
              clear auxFEVD
              %Note that we are computing the Cholesky FEVD for each draw,
              %eventhough we do not really need it. 
          end
      end
      
      clear AL vechSigma Gamma Csim Cauxsim Ccumsim B1 B1aux B1cholsim
   end
    display('The numerical procedure is over.')
    toc   

%clear ip;

%%

auxtxt1=100*mean(Inference.stat,2);

display(blanks(1));
display('Remarks on the numerical procedure:')
display(strcat('a)',num2str(auxtxt1),...
'% of the Draws have an stationary A'));

auxtxt2=100*mean(Inference.pdSigma,2);

display(blanks(1));
display(strcat('b)',num2str(auxtxt2),...
'% of the Draws have a positive definite Sigma'));

clear auxtxt1;


%% 9) Implement "Standard" Bootstrap Inference

   aux=reshape(Inference.stat.*Inference.pdSigma,[1,1,I]);
   %Inference.bootsIRFZ=quantile(Inference.IRFZ(:,:,aux==1,:),[.16,.84],3);
   Inference.bootsIRFZ=quantile(Inference.IRFZ(:,:,aux==1,:),[((1-confidence)/2),1-((1-confidence)/2)],3);
   Inference.bootsPointEstFEVD=quantile(Inference.PointEstFEVD(:,:,aux==1),[((1-confidence)/2),1-((1-confidence)/2)],3);
   %This would be a standard Efron bootstrap for FEVD 
   Inference.bootsFEVD=quantile(Inference.FEVD(:,:,aux==1,:),[((1-confidence)/2),1-((1-confidence)/2)],3);
   %This would give one sided confidence sets for the bounds on the FEVDs
   %described in the paper
%% 10) Plots of IRFs using external instruments
addpath('3AuxFunctions/figuresfun')

figure(1)
subplot(3,1,1)
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(1,:,2,2),Inference.bootsIRFZ(1,:,1,2),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h1=plot(Inference.IRFZ(1,:,I,2),'b'); hold on
plot(InferenceMSW.MSWuboundcum(1,:),'--b'); hold on
plot(InferenceMSW.MSWlboundcum(1,:),'--b'); hold on
plot(Chol(1,:,2),'r'); hold on
set(h1,'LineWidth',3); hold off
clear h1
axis([1 20 -40 40]); 
xlabel('Months after the shock');
ylabel('Percent');
title('Cumulative Response of the Growth rate of Oil Production'); 

subplot(3,1,2)
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(2,:,2,1),Inference.bootsIRFZ(2,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h2=plot(Inference.IRFZ(2,:,I,1),'b'); hold on
plot(InferenceMSW.MSWubound(2,:),'--b'); hold on
plot(InferenceMSW.MSWlbound(2,:),'--b'); hold on
plot(Chol(2,:,1),'r'); hold on
set(h2,'LineWidth',3); hold off
clear h2
axis([1 20 -40 40]);   
xlabel('Months after the shock');
ylabel('Percent');
%legend('95% "Bootstrap" Confidence Set (J=1,000)','Point Estimator Stock & Watson','95% MSW Confidence Set','Location','southwest');
title('Response of Global Real Activity');
clear h2;

subplot(3,1,3)
[~,~] = jbfill(1:1:hori+1,Inference.bootsIRFZ(3,:,2,1),Inference.bootsIRFZ(3,:,1,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h3=plot(Inference.IRFZ(3,:,I,1),'b'); hold on
plot(InferenceMSW.MSWubound(3,:),'--b'); hold on
plot(InferenceMSW.MSWlbound(3,:),'--b'); hold on
plot(Chol(3,:,1),'r'); hold on
set(h3,'LineWidth',3); hold off
clear h3
axis([1 20 -40 40]); 
xlabel('Months after the shock');
ylabel('Percent');
title('Response of the Real Price of Oil'); clear h3;

%% 11) Forecast-Error Variance Decompositions

figure(2)
subplot(3,1,1)
[~,~] = jbfill(1:1:hori+1,Inference.bootsPointEstFEVD(1,:,2),Inference.bootsPointEstFEVD(1,:,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h1=plot(Inference.PointEstFEVD(1,:,I),'b'); hold on
plot(Inference.CholFEVD(1,:,I),'r'); hold on
plot(Inference.bootsFEVD(1,:,1,2),'--b'); hold on
plot(Inference.bootsFEVD(1,:,2,1),'--b'); hold on
set(h1,'LineWidth',3); hold off
clear h1
axis([1 20 0 1]);  
xlabel('Months after the shock');
ylabel('Percent');
title('Contribution of the supply shock to the FE of the Growth rate of Oil Production');  

subplot(3,1,2)
[~,~] = jbfill(1:1:hori+1,Inference.bootsPointEstFEVD(2,:,2),Inference.bootsPointEstFEVD(2,:,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h2=plot(Inference.PointEstFEVD(2,:,I),'b'); hold on
plot(Inference.CholFEVD(2,:,I),'r'); hold on
plot(Inference.bootsFEVD(2,:,1,2),'--b'); hold on
plot(Inference.bootsFEVD(2,:,1,1),'--b'); hold on
set(h2,'LineWidth',3); hold off
clear h2;
axis([1 20 0 1]);  
xlabel('Months after the shock');
ylabel('Percent');
%legend('95% "Bootstrap" Confidence Set (J=1,000)','Point Estimator Stock & Watson','95% MSW Confidence Set','Location','southwest');
title('Contribution of the supply shock to the FE of Global Real Activity');

subplot(3,1,3)
[~,~] = jbfill(1:1:hori+1,Inference.bootsPointEstFEVD(3,:,2),Inference.bootsPointEstFEVD(3,:,1),[204/255 204/255 204/255],[204/255 204/255 204/255],0,0.5); hold on
h3=plot(Inference.PointEstFEVD(3,:,I),'b'); hold on
plot(Inference.CholFEVD(3,:,I),'r'); hold on
plot(Inference.bootsFEVD(3,:,1,2),'--b'); hold on
plot(Inference.bootsFEVD(3,:,1,1),'--b'); hold on
set(h3,'LineWidth',3); hold off
clear h3
axis([1 20 0 1]);  
xlabel('Months after the shock');
ylabel('Percent');
title('Contribution of the supply shock to the FE of the Real Price of Oil'); clear h3;



%% 12) Plot of the estimated target shock

%Create anualized value of the structural shocks:

[Yearaux,Monthaux] = meshgrid(1973:2004, 1:12);
Yearsaux=Yearaux(:); clear Yearaux;
Monthsaux=Monthaux(:); clear Monthaux;
Datesaux=[Yearsaux(26:end-3,1),Monthsaux(26:end-3,1),epsilonhat',SVARinp.Z(25:end,1)];
%Starts in 1975 because we have 2 years of lags

for yearx=1:30
shockannual(yearx,:)=sum(Datesaux(Datesaux(:,1)==1974+yearx,3))/sum(Datesaux(:,1)==1974+yearx);
instrumentannual(yearx,:)=sum(Datesaux(Datesaux(:,1)==1974+yearx,4))/sum(Datesaux(:,1)==1974+yearx);
end

figure(3)
for j=1:8
celldates(j,:)={dataU.dates(startdata+1+p+((j-1)*50),:)};
end

plot(1:T,epsilonhat,'b'); hold on
h1=plot([1 T],[std(epsilonhat) std(epsilonhat)],'black'); hold on
h2=plot([1 T],[2*std(epsilonhat) 2*std(epsilonhat)],'r'); hold on
plot([1 T],[-std(epsilonhat) -std(epsilonhat)],'black'); hold on
plot([1 T],[-2*std(epsilonhat) -2*std(epsilonhat)],'r'); hold off
xlim([1,T]);
legend([h1,h2],{'One standard deviation','Two standard deviations'})
ax=gca;
ax.XTick =1:50:351;
ax.XTickLabel=celldates;
ax.XTickLabelRotation = 90;
ylabel('Estimated Oil Shock');
grid on

figure(4)
plot(Datesaux(1,1):1:Datesaux(end,1),-shockannual)
grid on
axis([Datesaux(1,1) Datesaux(end,1) -1 1]);
xlabel('Year')



%% Save output   
cd('../6SVAR_MSW_pics');
namdir=strcat(date,'_',instrumentname,'_','p=',num2str(p),'confidence',num2str(confidence));
display(blanks(1));
display(strcat('Figures and output are in:',blanks(3),'SVAR_SW_pics\',namdir));
if exist(namdir,'dir')
    cd(namdir);
    print(figure(1),'-depsc2',strcat('IRFs_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(2),'-depsc2',strcat('FEVDs_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(3),'-depsc2',strcat('Shock_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(4),'-depsc2',strcat('Shock_SW_Annual_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    save(strcat('IRFs_SW_','p=',num2str(p),'_',num2str(startyear),'.mat'));
else
    mkdir(namdir)
    cd(namdir)
    print(figure(1),'-depsc2',strcat('IRFs_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(2),'-depsc2',strcat('FEVDs_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(3),'-depsc2',strcat('Shock_SW_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    print(figure(4),'-depsc2',strcat('Shock_SW_Annual_','p=',num2str(p),'_',num2str(startyear),'.eps'));
    save(strcat('IRFs_SW_','p=',num2str(p),'_',num2str(startyear),'.mat'));
end
cd ..
cd ..


