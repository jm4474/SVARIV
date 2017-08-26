function [c,der_c] = upperchol_id_minusj(AL,Sigma,Gamma)
%  -Zero restriction for the two shock-two instrument SVARIV model (-jth shock)
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

n                    = size(Sigma,1);

aux                  = eye(n);

Sigmainv             = (Sigma\aux);

[IRFSVARIV, dIRFdmu] = IRFSVARIV_2inst_j(AL,Sigma,Gamma,...
                      1,-1,2,@upperchol_id);

c                    = Sigmainv*IRFSVARIV(:,1);

%% 2) Define the vector of derivatives

p           = size(AL,2)/n;

dcdvecA     = dIRFdmu(1:(n^2)*p,:,1)*Sigmainv;

M           = vechtovec(n);

dcdvecSigma = dIRFdmu(((n^2)*p)+1:((n^2)*p)+(n*(n+1)/2),:,1)*Sigmainv... 
             -M'*kron(Sigmainv,Sigmainv)*kron(IRFSVARIV(:,1),aux);

dcdvecGamma = dIRFdmu(((n^2)*p)+(n*(n+1)/2)+1:end,:,1)*Sigmainv;

der_c        = [ dcdvecA; dcdvecSigma; dcdvecGamma ];

end

