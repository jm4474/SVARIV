function [ ARteststatistic, Gamma11 ] = ARteststatistic(AL,Gamma,horizons,scale,norm)
%  -Computes IRFs identified using an external instrument
%  -Syntax:
%    [ ARteststatistic ] = ARteststatistic(AL,Gamma,horizons,scale,norm, lambda)
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

n         = size(Gamma,1);

p         = size(AL,2)/n;

%Reduced-form MA coefficients
Cauxsim   = [eye(n),MARep(AL,p,horizons)]; 
      
Csim      = reshape(Cauxsim,[n,n,horizons+1]);

%Ccumsim   = cumsum(Csim,3);

Gamma11 = Gamma(norm,1);

ARteststatistic(:,:) = reshape(sum(bsxfun(@times,Csim,(Gamma*scale)'),2),[n,horizons+1]);


%ARteststatistic(:,:,1) = reshape(sum(bsxfun(@times,Csim,(Gamma*scale)'),2),[n,horizons+1]);
%ARteststatistic(:,:,2) = reshape(sum(bsxfun(@times,Ccumsim,(Gamma*scale)'),2),[n,horizons+1]);



      