function [mu,AL,Sigma,eta,X,Y] = RForm_VAR(TSL,p,W)
% -------------------------------------------------------------------------
% Provides reduced form estimators of a VAR(p) model  
% 
% Inputs:
% - TSL: matrix of dimension (T times n) containing the time series
% - p: number of lags in the VAR model
% Outputs:
% - AL: VAR model coefficients
% - Sigma: covariance matrix of VAR model residuals
% - eta: VAR model residuals
% - X: VAR model regressors
% - W: Matrix of dummies
%
% The estimation always includes a constant
% 
% This version: May 4th, 2015
% Last edited by Jos?Luis Montiel-Olea

% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------


%% Definitions
aux = lagmatrix(TSL,1:1:p);
Y= TSL((p+1):end,:); %The rows of this matrix are Y_t'
X= [ones(size(Y,1),1),aux((p+1):end,:)]; clear aux   %The rows of this matrix are [1,X_{t}'] in p.6
if nargin==3
slopeparameters=(Y'*W*X)*((X'*W*X)^(-1)); %contains the nx1 constant vector and AL
else
slopeparameters=(Y'*X)*((X'*X)^(-1));
end
%%%%%%%%%%COMMENT: MAY 4th
%%%%%%%%%%The model now includes a constant. 

%% Generate \mu and the vec(A(L)) estimators
AL=slopeparameters(:,2:end); %n x np
mu=slopeparameters(:,1);  %n x 1


%% Covariance matrix
eta = Y'-slopeparameters*(X'); 
Sigma = (eta*eta')/(size(eta,2));


end

