function [IRFs, bootsIRFs] = Gasydistboots(seed, I, n, p, nvar, x, hori, confidence, T, vecAL, vechSigma, Gamma, Whatall, f)
%  -Provides inference for SVAR-IV based on samples from the asy. dist.
%  -Syntax:
%    [IRFs, bootsIRFs] = Gasydistboots(seed, I, n, p, nvar, x, hori, confidence, vecAL, vechSigma, Gamma, Whatall)
%  -Inputs:
%     seed: seed structure  
%        I: number of samples from the asymptotic distribution
%        n: number of variables in the VAR
%        p: number of lags in the VAR
%     nvar: normalizing variable
%        x: scale
%     hori: horizons
%confidence: confidence level
%        T: time periods
%        f: function handle (depends on AL, Sigma, Gamma, hori, x, nvar)
%    vecAL: point estimator of vec(AL)
%vechSigma: point estimator of vech(Sigma)
%    Gamma: point estimator of vec(Gamma)
%  Whatall: Covariance mat of (vec(A)', vech(Sigma)', vec(Gamma)')              
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

%% 1) Set the seed

rng(seed); clear seed

%% 2) Make sure that Whatall is symmetric and positive semidefinite

dall        = size(Whatall,1);

Whatall     = (Whatall + Whatall')/2;
    
[aux1,aux2] = eig(Whatall);
    
Whatall     = aux1*max(aux2,0)*aux1'; 

%% 3) Generate draws from the Gaussian asy. dist. 
% (centered at point estimators)

gvar    = [mvnrnd(zeros(I,dall),(Whatall)/T)',...
                     zeros(dall,1)];
          %Added an extra column of zeros to access point estimators       
    
Draws   = bsxfun(@plus,gvar,...
          [vecAL;vechSigma;Gamma(:)]);

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

IRFs       = zeros(n, hori+1, ndraws);
    
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
          
    IRFs(:,:,idraws) = f(AL,Sigma,Gamma,hori,x,nvar);
            
    clear AL vechSigma Gamma 
      
end
   
%% 5) Implement "Standard" Bootstrap Inference

aux        = reshape(pdSigma,[1,1,ndraws]);
   
bootsIRFs  = quantile(IRFs(:,:,aux==1),...
                          [((1-confidence)/2),1-((1-confidence)/2)],3);      

end
