function [IRFs, bootsIRFs] = Gasydistboots(seed, NB, n, p, norm, scale, horizons, confidence, T, f, SVARinp, NWlags, AL, Sigma, Gamma, V, WHatall)
%  -Provides inference for SVAR-IV based on samples from the asy. dist.
%  -Syntax:
%    [IRFs, bootsIRFs] = Gasydistboots(seed, NB, n, p, norm, scale, horizons, confidence, T, f, SVARinp, NWlags, AL, Sigma, Gamma, V, Whatall)
%  -Inputs:
%     seed: seed structure  
%       NB: number of samples from the asymptotic distribution
%        n: number of variables in the VAR
%        p: number of lags in the VAR
%     norm: normalizing variable
%    scale: scale of the shock
% horizons: number of horizons (IRFs) (does not include the impact or horizon 0)
%confidence: confidence level
%        T: time periods
%        f: function handle (depends on AL, Sigma, Gamma, hori, x, nvar)
%   NWlags: Newey-West lags
%       AL: point estimator of AL
%vechSigma: point estimator of vech(Sigma)
%    Gamma: point estimator of vec(Gamma)
%  Whatall: Covariance mat of (vec(AL)', vech(Sigma)', vec(Gamma)')              
%
%  -Output:
%       AL: Least-squares estimator of the VAR coefficients
%    Sigma: Least-squares estimator of the VAR residuals
%      eta: VAR model residuals
%        X: VAR model regressors
%        Y: VAR matrix of endogenous regressors
%  -Output:
%     IRFs: 3d structure containing all of the bootstrap draws of IRF
% bootsIRF: alpha/2 and 1-alpha/2 quantiles of IRFs
% 
% This version: July 17th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Create an RForm (if necessary)


%a) Estimation of (AL, Sigma) and the reduced-form innovations
if (nargin == 12)

    % This essentially checks whether the user provided or not his RForm. If
    % the user didn't, then we calculate it. If the user did, we skip this section and
    % use his/her RForm.
    
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
        
    [RForm.WHatall,RForm.WHat,RForm.V] = ...
    CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,NWlags);                
    

elseif (nargin == 17)

    % We will use the user's RForm

else 

    % Error, not enough inputs.

end

vechSigma = V*Sigma(:);

AL = AL(:);

%% 2) Set the seed

rng(seed); clear seed

%% 2) Make sure that Whatall is symmetric and positive semidefinite

dall        = size(WHatall,1);

WHatall     = (WHatall + WHatall')/2;
    
[aux1,aux2] = eig(WHatall);
    
WHatall     = aux1*max(aux2,0)*aux1'; 

%% 3) Generate draws from the Gaussian asy. dist. 
% (centered at point estimators)

gvar    = [mvnrnd(zeros(NB,dall),(WHatall)/T)',...
                     zeros(dall,1)];
          %Added an extra column of zeros to access point estimators       
    
Draws   = bsxfun(@plus,gvar,...
          [AL;vechSigma;Gamma(:)]);

k       = size(Gamma,1)/n;      

%The vector "Draws" represents a vector of I draws
%from a multivariate normal vector (of dimension dall) centered
%at [vec(A)', vech(Sigma)',Gamma']' with covariance matrix 
%(WHatall/T). Thus, it represents a draw from the asy. dist
%of the reduced-form parameters.      

%% 4) Evaluate the parameter of interest 
% (which is allowed to depend on the full vector vecA, vechSigma, Gamma)

ndraws     = size(Draws,2);
     
pdSigma    = zeros(1,ndraws);

IRFs       = zeros(n, horizons+1, ndraws, 2);

for idraws = 1:ndraws
    
      %i) Generate the draws for AL 
      
      AL   = reshape(Draws(1:(n^2)*p,idraws),[n,n*p]);
      
      %ii) Generate the draws from Sigma
      
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
          
    [IRFs(:,:,idraws,:),~] = f(AL,Sigma,Gamma,horizons,scale,norm);
            
    clear AL vechSigma Gamma 
      
end
   
%% 5) Implement "Standard" Bootstrap Inference

aux        = reshape(pdSigma,[1,1,ndraws]);
   
bootsIRFs  = quantile(IRFs(:,:,aux==1,:),...
                          [((1-confidence)/2),1-((1-confidence)/2)],3);      

end
