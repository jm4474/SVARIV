function [ IRFSVARIV, Gamma11 ] = IRFSVARIV(AL,Sigma,Gamma,horizons,scale,norm)
%  -Computes IRFs identified using an external instrument
%  -Syntax:
%    [ IRFSVARIV ] = IRFSVARIV(AL,Gamma,hori,x,nvar)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Gamma: covariances between zt and etat              (n x 1)
%     horizons: number of horizons to evaluate IRFs      (1 x 1)  
%        scale: scale of the shock                       (1 x 1)
%     norm: normalizing variable                         (1 x 1)
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

%Reduced-form MA coefficients
Cauxsim   = [eye(n),MARep(AL,p,horizons)]; 
      
Csim      = reshape(Cauxsim,[n,n,horizons+1]);

Ccumsim   = cumsum(Csim,3);
      
B1        = scale*Gamma./Gamma(norm,1);    
      
IRFSVARIV(:,:,1) = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,horizons+1]);
IRFSVARIV(:,:,2) = reshape(sum(bsxfun(@times,Ccumsim,B1'),2),[n,horizons+1]);
%Ck(a)*x*Gamma/Gamma11. Third dimension for cumulative and non-cumulative


Gamma11 = Gamma(norm,1);

end


