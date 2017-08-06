function [c,der_c] = lowerchol_id(AL,Sigma,Gamma)
%  -Zero restriction for the two shock-two instrument SVARIV model
%  -Syntax:
%    [ c,der_c,e_j ] = zerorestriction(A,Sigma,Gamma)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Sigma: covariance matrix of residuals               (n x n)
%    Gamma: covariance matrix between inst and shocks    (n x 2)
%  -Output:
%        c: vector defining the zero restriction         (n x 1)
%    der_c: matrix of derivatives of c                   (n^2p + (n*(n+1)/2) + 2n x n)
%
% This version: August 4th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Define the vector of zero restrictions implied by the lower Chol decomposition

n           = size(Sigma,1);

aux         = eye(n);

Sigmainv    = (Sigma\aux);

c           = Sigmainv*aux(:,2);

%% 2) Define the vector of derivatives

p           = size(AL,2)/n;

sizeGamma   = size(Gamma(:),1);

dcdvecA     = zeros((n^2)*p,n);

M           = vechtovec(n);

dcdvecSigma = -M'*kron(Sigmainv,Sigmainv)*kron(aux(:,2),aux);

dcdvecGamma = zeros(sizeGamma,n);

der_c        = [ dcdvecA; dcdvecSigma; dcdvecGamma ];


end

