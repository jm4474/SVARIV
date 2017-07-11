function [mu,AL,Sigma,eta,X,Y] = RForm_VAR(TSL,p,W)
% -Provides reduced form estimators of a VAR(p) model 
% -Syntax:
%       [mu,AL,Sigma,eta,X,Y] = RForm_VAR(TSL,p,W)
% -Inputs:
%     TSL: matrix of time series                       (T times n)
%       p: number of lags in the VAR model             (1 times 1)
%       W: weighting matrix (possibly for dummies)     (T times T)
% -Output:
%      AL: Least-squares estimator of the VAR coefficients
%   Sigma: Least-squares estimator of the VAR residuals
%     eta: VAR model residuals
%       X: VAR model regressors
%       Y: VAR matrix of endogenous regressors
%   
% -Note  : estimation always includes a constant
% 
% This version: July 11th, 2017
% Last edited by José Luis Montiel-Olea

%% Definitions

aux  = lagmatrix(TSL,1:1:p);

Y    = TSL((p+1):end,:); 
%The rows of this matrix are Y_t'

X    = [ones(size(Y,1),1),aux((p+1):end,:)]; clear aux   
%The rows of this matrix are [1,X_{t}'] 

if nargin == 3
    
    slopeparameters = (Y'*W*X)*((X'*W*X)^(-1)); 

else
    
    slopeparameters = (Y'*X)*((X'*X)^(-1));
end

%% Generate \mu and vec(A(L)) estimators

AL   = slopeparameters(:,2:end); %n x np

mu   = slopeparameters(:,1);  %n x 1

%% Covariance matrix

eta  = Y'-slopeparameters*(X'); 

Sigma= (eta*eta')/(size(eta,2));

end
