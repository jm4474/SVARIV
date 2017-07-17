function [ IRFSVARIV ] = IRFSVARIV(AL,Sigma,Gamma,hori,x,nvar)
%  -Computes IRFs identified using an external instrument
%  -Syntax:
%    [ IRFSVARIV ] = IRFSVARIV(AL,Gamma,hori,x,nvar)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Gamma: covariances between zt and etat              (n x 1)
%     hori: number of horizons to evaluate IRFs          (1 x 1)  
%        x: scale                                        (1 x 1)
%     nvar: normalizing variable                         (1 x 1)
%  -Output:
%IRFSVARIV: vector of IRFs                               (n x hori+1)    
%   
% -Note  : estimation always includes a constant
% 
% This version: July 17th, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Main definitions and function

n         = size(Sigma,1);

p         = size(AL,2)/n;

Cauxsim   = [eye(n),MARep(AL,p,hori)]; 
      
Csim      = reshape(Cauxsim,[n,n,hori+1]);
      
B1        = x*Gamma./Gamma(nvar,1);    
      
IRFSVARIV = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);

end

