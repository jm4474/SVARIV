function [c,der_c] = lowerchol_id(AL,Sigma)
%  -Zero restriction for the two shock-two instrument SVARIV model
%  -Syntax:
%    [ c,der_c,e_j ] = zerorestriction(A,Sigma,j)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Sigma: covariance matrix of residuals               (n x n)
%        j: index of the restricted column               (1 x 1)
%           (j is either 1 or 2)
%  -Output:
%        c: vector defining the zero restriction         (n x 1)
%    der_c: vector of derivatives of c                   (n x 1)
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

dcdvecA     = zeros((n^2)*p,n);

M           = vechtovec(n);

dcdvecSigma = -M'*kron(Sigmainv,Sigmainv)*kron(aux(:,2),aux);

der_c        = [ dcdvecA; dcdvecSigma ];


end

