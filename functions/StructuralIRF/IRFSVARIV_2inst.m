function [ IRFSVARIV ] = IRFSVARIV_2inst(AL,Sigma,Gamma,hori,x,nvar)
%  -Computes IRFs identified using two external instruments possibly
%  correlated with the two target shocks
%  -Syntax:
%    [ IRFSVARIV ] = IRFSVARIV(AL,Gamma,hori,x,nvar)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Gamma: covariances between zt and etat              (n x 2)
%     hori: number of horizons to evaluate IRFs          (1 x 1)  
%        x: scale                                        (1 x 1)
%     nvar: normalizing variable                         (1 x 1)
%  -Output:
%IRFSVARIV: vector of IRFs                               (n x hori+1)    
%   
% 
% This version: July 17th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Main definitions and function

n         = size(Sigma,1);

p         = size(AL,2)/n;

Cauxsim   = [eye(n),MARep(AL,p,hori)]; 
      
Csim      = reshape(Cauxsim,[n,n,hori+1]);

Gamma1    = Gamma(:,1);

Gamma2    = Gamma(:,2);

aux1      = Sigma\Gamma1;

aux2      = Sigma\Gamma2;

Gamma3    = (Gamma1-((aux1(2,1)/aux2(2,1))*Gamma2));
      
B1        = x*Gamma3./Gamma3(nvar,1);    
      
IRFSVARIV = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);

end

