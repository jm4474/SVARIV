    function [WHataux,WHat,V] = CovAhat_Sigmahat_Gamma(p,X,Z,eta,lags)
% -------------------------------------------------------------------------
% Computes the asymptotic variance of [vec(Ahat)',vech(Sigmahat)',vec(Gammahat)']'
% 
% Inputs:
% - p: lag order
% - Sigma: covariance matrix of residuals
% - X: T times np + 1 (n is the dimension of the VAR)
% - Z: T times k      (k is the number of instruments)
% - eta: reduced-form residuals
% Outputs:
% - WHat: asymptotic variance of [vec(Ahat)',vech(Sigmahat)',vec(Gammahat)']'
% (this matrix also includes Sigmahat, altough we do not use this info) 
%
%
% This version: October 6th, 2016
% Last revised by Jos�-Luis Montiel Olea


% -------------------------------------------------------------------------


%% Definitions
n = size(eta,1);
k = size(Z,2);
XSVARp = X; 


matagg = [XSVARp,eta',Z]'; %The columns of this vector are (1;X_t; eta_t;Z_t)
T1aux = size(eta,2); %This is the number of time periods
T2aux = size(matagg,1); %This is the column dimension of (1;X_t;eta_t;Z_t)
etaaux = reshape(eta,[n,1,T1aux]); %Each 2-D page contains eta_t
mataggaux = permute(reshape(matagg,[T2aux,1,T1aux]),[2,1,3]); %Each 2-D page contains (1,X_t',\eta_t',Z_t')
auxeta = bsxfun(@plus,bsxfun(@times,etaaux,mataggaux),-mean(bsxfun(@times,etaaux,mataggaux),3));
%Each 2-D page contains [eta_t, eta_t X_t', eta_t eta_t'-Sigma, eta_tZ_t'-Gamma];
vecAss1= reshape(auxeta,[n+(p*(n^2))+n^2+(n*k),1,T1aux]);
%Each 2-D page contains [eta_t; vec(eta_tX_t') ; vec(eta_t*eta_t'-Sigma) ; vec(eta_tZ_t'-Gamma)]
WhatAss1 = sum(bsxfun(@times,vecAss1,permute(vecAss1,[2,1,3])),3)./T1aux;
%This is the covariance matrix we are interested in 
%Zaux = reshape(Z,[k,1,T1aux]);
%NWhatAss = (1-1/(8+1))*sum(Zaux(:,:,1:end-1).*Zaux(:,:,2:end).*...
        %bsxfun(@times,etaaux(:,:,1:end-1),permute(etaaux(:,:,2:end),[2,1,3])),3);
%WhatAss2 = WhatAss1
for i_n = 1:lags
    %NWhatAss = sum(Zaux(:,:,1:end-i_n).*Zaux(:,:,i_n+1:end).*...
        %bsxfun(@times,etaaux(:,:,1:end-i_n),permute(etaaux(:,:,i_n+1:end),[2,1,3])),3);
    NWhatAss = (1/T1aux-i_n)*...
        sum(bsxfun(@times,vecAss1(253:261,:,1:end-i_n),permute(vecAss1(253:261,:,i_n+1:end),[2,1,3])),3);
    for i = 1:9
        for j = 1:9
            WhatAss1(252+i,252+j) = WhatAss1(252+i,252+j) + (1-i_n/(lags+1))*(NWhatAss(i,j)+NWhatAss(j,i));
        end
    end
    clear NWhatAss 
end
clear i j i_n

%% Construct the selector matrix Vaux that gives: vech(Sigma)=Vaux*vec(Sigma)
I = eye(n);
V = kron(I(1,:),I);
for i=2:n
    V = [V; kron(I(i,:),I(i:end,:))];
end


%% This is the estimator for matrix What MSW (2015)
Q1=(XSVARp'*XSVARp./T1aux);
Q2=Z'*XSVARp/T1aux;
Shat = [kron([zeros(n*p,1),eye(n*p)]*(Q1^(-1)),eye(n)), zeros((n^2)*p,n^2+(k*n)); ...
        zeros(n*(n+1)/2,((n^2)*p)+n), V, zeros(n*(n+1)/2,k*n);...
         -kron(Q2*(Q1^(-1)),eye(n)), zeros(k*n,n^2),eye(k*n) ];      
WHataux = (Shat)*(WhatAss1)*(Shat');

%WHataux is the covariance matrix of vec(A),vech(Sigma),Gamma

WHat=[WHataux(1:(n^2)*p,1:(n^2)*p),...
    WHataux(1:(n^2)*p,((n^2)*p)+(n*(n+1)/2)+1:end);...
    WHataux(1:(n^2)*p,((n^2)*p)+(n*(n+1)/2)+1:end)',...
    WHataux(((n^2)*p)+(n*(n+1)/2)+1:end,((n^2)*p)+(n*(n+1)/2)+1:end)];

%WHat is the part of WHataux corresponding to vec(A) and Gamma
    
%%%%%%%%%%%%%%COMMENT 
%%%%%%%%%%%%%%Note that when the number of parameters in the covariance
%%%%%%%%%%%%%%matrix is large, relative to n, the matrix will not be
%%%%%%%%%%%%%%invertible. In particular, we need to guarantee that
%%%%%%%%%%%%%% n^2 p + n(n+1)/2 + nK is strictly smaller than T. Otherwise, the
%%%%%%%%%%%%%% covariance matrix is not invertible. We should display this as
%%%%%%%%%%%%%% a warning message. Even when the matrix is invertible, we
%%%%%%%%%%%%%% have to be careful with the conditioning number of the
%%%%%%%%%%%%%% matrix

end
