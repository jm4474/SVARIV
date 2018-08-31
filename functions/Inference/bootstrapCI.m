%% 1) Set number of VAR Lags, Newey West lags and confidence level.

p           = 24;     %Number of lags in the VAR model
 
confidence  = 1;    %Confidence Level for the standard and weak-IV robust confidence set

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
                 
NB          = 1000; % number of samples from the asymptotic distribution

%% 2) Load data (saved in structure "data")
%  These are the variables that were defined in line 14 above. 
%  The time units should be on a single.xls file called Time.xls
%  All the VAR variables should be on a single .xls file called Data.xls
%  The external instrument should be in a single .xls file called ExtIV.xls

cd ..

cd ..

direct = pwd;

application = 'Oil';

cd(strcat(direct,'/Data/',application));

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


%% 3) Create an RForm (if necessary)

SVARinp.ydata = ydata;

SVARinp.Z = z;

SVARinp.n        = size(ydata,2); %number of columns(variables)

RForm.p          = p; %RForm_user.p is the number of lags in the model


%a) Estimation of (AL, Sigma) and the reduced-form innovations


% This essentially checks whether the user provided or not his RForm. If
% the user didn't, then we calculate it. If the user did, we skip this section and
% use his/her RForm.

addpath('functions/RForm');

[RForm.mu, ...
 RForm.AL, ...
 RForm.Sigma,...
 RForm.eta,...
 RForm.X,...
 RForm.Y]        = RForm_VAR(SVARinp.ydata, p);

%RForm.AL(:),RForm.V*RForm.Sigma(:),RForm.Gamma(:), RForm.WHatall

%b) Estimation of Gammahat (n times 1)

RForm.Gamma      = RForm.eta*SVARinp.Z(p+1:end,1)/(size(RForm.eta,2));   %sum(u*z)/T. Used for the computation of impulse response.
%(We need to take the instrument starting at period (p+1), because
%we there are no reduced-form errors for the first p entries of Y.)

%c) Add initial conditions and the external IV to the RForm structure

RForm.Y0         = SVARinp.ydata(1:p,:);

RForm.externalIV = SVARinp.Z(p+1:end,1);

RForm.n          = SVARinp.n;

%a) Covariance matrix for vec(A,Gammahat). Used
%to conduct frequentist inference about the IRFs. 
[RForm.WHatall,RForm.WHat,RForm.V] = ...
CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,NWlags);     

%NOTES:
%The matrix RForm.WHatall is the covariance matrix of 
% vec(Ahat)',vech(Sigmahat)',Gamma')'
 
%The matrix RForm.WHat is the covariance matrix of only
% vec(Ahat)',Gamma')' 
 
% The latter is all we need to conduct inference about the IRFs,
% but the former is needed to conduct inference about FEVDs.

vechSigma = RForm.V * RForm.Sigma(:);

AL = RForm.AL(:);

Gamma = RForm.Gamma;

%% 4) Estimation of the asymptotic variance of A,Gamma

% Definitions

n            = RForm.n; % Number of endogenous variables

T            = (size(RForm.eta,2)); % Number of observations (time periods)

d            = ((n^2)*p)+(n);     %This is the size of (vec(A)',Gamma')'

%dall         = d+ (n*(n+1))/2;    %This is the size of (vec(A)',vec(Sigma), Gamma')'



%% 5) Make sure that Whatall is symmetric and positive semidefinite

%dall        = size(RForm.WHatall,1);  


dall        = size(RForm.WHatall,1);

WHatall     = (RForm.WHatall + RForm.WHatall')/2;
    
[aux1,aux2] = eig(WHatall);
    
WHatall     = aux1*max(aux2,0)*aux1'; 

%% 6) Generate draws
% Centered at (vec(AL)', Gamma')'


gvar    = [mvnrnd(zeros(NB,dall),(WHatall)/T)',...
                     zeros(dall,1)];
          %Added an extra column of zeros to access point estimators       
    
Draws   = bsxfun(@plus,gvar,...
          [AL;vechSigma;Gamma(:)]);

k       = size(Gamma,1)/n;        

%% 7 Evaluate the parameter of interest  

f = @ARteststatistic;

ndraws     = size(Draws,2);
     
pdSigma    = zeros(1,ndraws);

IRFs       = zeros(n, horizons+1, ndraws);

Gamma11s   = zeros(1,ndraws);

addpath('functions/StructuralIRF');

RFormIRFBoots = zeros(n, horizons + 1,ndraws);

AlphaBoots = zeros(1, ndraws);

for idraws = 1:ndraws
        
        %i) Generate the draws for AL 
        
        AL   = reshape(Draws(1:(n^2)*p,idraws),[n,n*p]);
        
        %ii) Draws from Sigma
      
        vechSigma = Draws((n^2)*p+1:(n^2)*p+(n*(n+1)/2),idraws);
      
        Sigma = tril(ones(n),0); Sigma(Sigma==1) = vechSigma';
      
        Sigma = Sigma + tril(Sigma,-1)';
      
        %Check if the draws are positive definite
      
          if min(eig(Sigma))>0
          
          pdSigma(1,idraws) = 1;
          
          else
          
          pdSigma(1,idraws) = 0;
          
          end
          
          %iii) Draws from Gamma
      
          Gamma = reshape(Draws(((n^2)*p)+(n*(n+1)/2)+1:end,idraws),...
              [n,k]);  
          
          [RFormIRFBoots(:,:,idraws), AlphaBoots(idraws)] = f(AL,Gamma,horizons,scale,norm);
            
      clear AL vechSigma Gamma 
    
end

grid = rand(100,1);

grid_size = size(grid,1);

difference       = zeros(grid_size, ndraws, n, horizons+1);

for var = 1:n

    for horizon = 1:horizons+1
        
        IRFBootsVH = RFormIRFBoots(var,horizon,:);
        
        IRFBootsVH = reshape(IRFBootsVH, 1, 1001);
          
        difference(:,:,var,horizon) = (IRFBootsVH - (grid * AlphaBoots(1,:)))*(T^.5);
        %Eventually this has to be the AR statistic
        
        clear IRFBootsVH;
        
    end
    
end

%THIS IS A HORRIBLE NAME. WE WILL CHANGE IT. 
new_difference = (difference - difference(:,1001,:,:)).^2;
%Adjust so that the number of draws can vary. the element I pick should be
%ndraws

%% 8) Implement "Standard" Bootstrap Inference

aux        = reshape(pdSigma,[1,ndraws]);


bootsIRFs  = quantile(new_difference(:,aux==1,:,:),...
                          [((1-confidence)/2),1-((1-confidence)/2)],2);      
                      
difference_T = difference(:,1001,:,:);

% check whether each one would be rejected or not

% We are using a loop here, but the most efficient way would be to use
% logical indexing

%for var = 1:n
    
%    for horizon = 1:horizons+1
        
%        for lbd = 1:lambda
           
%            if(newest(lbd,1,var,horizon) >= bootsIRFs(lbd,1,var,horizon) && newest(lbd,1,var,horizon) <= bootsIRFs(lbd,2,var,horizon))
               
%                reject(n,horizon,lbd) = 0;
                
%            else
                
%                reject(n,horizon,lbd) = 1;
           
%            end
            
%        end
        
%    end
    
%end

reject = (difference_T(:,1,:,:) < bootsIRFs(:,1,:,:)) | (difference_T(:,1,:,:) > bootsIRFs(:,2,:,:));

reject = reshape(reject,[size(grid,1), n, horizons+1]);

