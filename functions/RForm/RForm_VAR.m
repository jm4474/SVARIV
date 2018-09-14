function [mu,AL,Sigma,eta,X,Y] = RForm_VAR(TSL,p,W)
% -Provides reduced form estimators of a VAR(p) model 
% -Syntax:
%       [mu,AL,Sigma,eta,X,Y] = RForm_VAR(TSL,p,W)
% -Inputs:
%     TSL: matrix of time series                           (T x n)
%       p: number of lags in the VAR model                 (1 x 1)
%       W: Matrix of exogenous regressors                  (T x m)
%          (could be dummies, deterministic time trends, etc)
% -Output:
%      AL: Least-squares estimator of the VAR coefficients (n x np)
%   Sigma: Least-squares estimator of the VAR residuals    (n x n)
%     eta: VAR model residuals                             (n x T)
%       X: Matrix of VAR covariates
%       Y: VAR matrix of endogenous regressors
%   
% -Note  : If no exogenous regressors are specified, our estimation 
%          always includes a constant.
% 
% This version: July 29th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) The lags of Y_t and Y_t itself

aux  = lagmatrix(TSL,1:1:p);

Y    = TSL((p+1):end,:); 
%The rows of this matrix are Y_t'

%% 2) Matrix of VAR Covariates
if nargin == 3

    X    = [W(p+1:end,:),aux((p+1):end,:)]; clear aux  
    
    m    = size(W,2);
%The rows of this matrix are [W_{t}',X_{t}'] 

else
    
    X    = [ones(size(Y,1),1),aux((p+1):end,:)]; clear aux  
    
    m    = 1;
%The rows of this matrix are [1,X_{t}']     
    
end    

%% 3) Estimate the coefficients for all the VAR covariates 

slopeparameters = (Y'*X)*((X'*X)^(-1));

% Alternative backslash implementation:

% slopeparameters = (X'*X\X'*Y)'

%% Generate the vec(A(L)) estimators 

AL   = slopeparameters(:,m+1:end); %n x np

mu   = slopeparameters(:,1:m);     %n x m (coefficients on the ex. reg.)

%% Covariance matrix

eta  = Y'-slopeparameters*(X'); 

Sigma= (eta*eta')/(size(eta,2));

end
