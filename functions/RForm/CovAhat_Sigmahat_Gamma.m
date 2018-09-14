function [WHataux,WHat,V] = CovAhat_Sigmahat_Gamma(p,X,Z,eta,lags)
% -Estimates the asymptotic covariance matrix of vec(A), vech(Sigma), Gamma
% -Syntax:
%       [WHataux,WHat,V] = CovAhat_Sigmahat_Gamma(p,X,Z,eta,lags)
% -Inputs:
%       p: VAR lags                                    (1 x 1)
%       X: VAR "right-hand" variables                  (T x (np+w))
%          (w is the number of additional regs)
%       Z: external instrument                         (T x k)
%     eta: eta                                         (n x T)
%    lags: Newey-West lags                             (1 x 1)
% -Output:
% WHataux: asymptotic variance of [vec(Ahat)',vech(Sigmahat)',vec(Gammahat)']'
%    WHat: asymptotic variance of [vec(Ahat)',vec(Gammahat)']'
%       V: Matrix such that vech(Sigma)=Vvec(Sigma) 
% 
% This version: July 29th, 2017
% Last revised by José-Luis Montiel Olea

%% Definitions
n      = size(eta,1);

k      = size(Z,2);

m      = size(X,2) - (n*p);

XSVARp = X; 


matagg = [XSVARp,eta',Z]'; %The columns of this vector are (W_t;X_t; eta_t;Z_t)

T1aux  = size(eta,2); %This is the number of time periods

T2aux  = size(matagg,1); %This is the column dimension of (W_t;X_t;eta_t;Z_t)

etaaux = reshape(eta,[n,1,T1aux]); %Each 2-D page contains eta_t

mataggaux = permute(reshape(matagg,[T2aux,1,T1aux]),[2,1,3]); %Each 2-D page contains (W_t',X_t',\eta_t',Z_t')

auxeta = bsxfun(@plus,bsxfun(@times,etaaux,mataggaux),-mean(bsxfun(@times,etaaux,mataggaux),3));

%Each 2-D page contains [eta_t W_t', eta_t X_t', eta_t eta_t'-Sigma, eta_tZ_t'-Gamma];

vecAss1 = reshape(auxeta,[(n*m)+(p*(n^2))+n^2+(n*k),1,T1aux]);

%Each 2-D page contains [vec(eta_tW_t'); vec(eta_tX_t') ; vec(eta_t*eta_t'-Sigma) ; vec(eta_tZ_t'-Gamma)]

% Auxiliary matrix to compute the HAC covariance matrix

AuxHAC1  = vecAss1(1:end,:,:);

AuxHAC2  = reshape(AuxHAC1,[size(vecAss1,1),size(vecAss1,3)])';

AuxHAC3  = NW_hac_STATA(AuxHAC2,lags);

WhatAss1 = AuxHAC3; 


%% Construct the selector matrix Vaux that gives: vech(Sigma)=Vaux*vec(Sigma)

I        = eye(n);

V        = kron(I(1,:),I);

for     i_vars=2:n
    V    = [V; kron(I(i_vars,:),I(i_vars:end,:))];
end


%% This is the estimator of What based on Montiel-Olea, Stock, and Watson

Q1       = (XSVARp'*XSVARp./T1aux);

Q2       = Z'*XSVARp/T1aux;

Shat     = [kron([zeros(n*p,m),eye(n*p)]*(Q1^(-1)),eye(n)), zeros((n^2)*p,n^2+(k*n)); ...
        zeros(n*(n+1)/2,((n^2)*p)+(n*m)), V, zeros(n*(n+1)/2,k*n);...
         -kron(Q2*(Q1^(-1)),eye(n)), zeros(k*n,n^2),eye(k*n) ];  
     
WHataux  = (Shat)*(WhatAss1)*(Shat');

%WHataux is the covariance matrix of vec(A),vech(Sigma),Gamma

WHat     = [WHataux(1:(n^2)*p,1:(n^2)*p),...
    WHataux(1:(n^2)*p,((n^2)*p)+(n*(n+1)/2)+1:end);...
    WHataux(1:(n^2)*p,((n^2)*p)+(n*(n+1)/2)+1:end)',...
    WHataux(((n^2)*p)+(n*(n+1)/2)+1:end,((n^2)*p)+(n*(n+1)/2)+1:end)];

%WHat is the part of WHataux corresponding to vec(A) and Gamma
    
%COMMENT 
%Note that when the number of parameters in the covariance
%matrix is large, relative to n, the matrix will not be
%invertible. In particular, we need to guarantee that
%n^2 p + n(n+1)/2 + nK is strictly smaller than T. 

end
