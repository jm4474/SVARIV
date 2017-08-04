function [ IRFSVARIV ] = IRFSVARIV_2inst_karel_pepesalgebra(AL,Sigma,Gamma,hori,x,nvar)
%  -Computes IRFs identified using two external instruments possibly
%  correlated with the two target shocks
%  -Syntax:
%    [ IRFSVARIV ] = IRFSVARIV(AL,Gamma,hori,x,nvar)
%  -Inputs:     
%       AL: matrix of autoregressive coefficients        (n x np)
%    Gamma: covariances between zt and etat              (n x k)
%     hori: number of horizons to evaluate IRFs          (1 x 1)  
%        x: scale                                        (1 x 1)
%     nvar: normalizing variable                         (1 x 1)
%  -Output:
%IRFSVARIV: vector of IRFs                               (n x hori+1)    
%   
% 
% This version: August 1st, 2017
% Last edited by José Luis Montiel-Olea

%% 1) Main definitions and function

n         = size(Sigma,1);

p         = size(AL,2)/n;

k         = size(Gamma,2);

Cauxsim   = [eye(n),MARep(AL,p,hori)]; 
      
Csim      = reshape(Cauxsim,[n,n,hori+1]);

%% 2) Partition Gamma according to Karel's AER notation

Sigmamu1 = Gamma(1:k,1:k);

Sigmamu2 = Gamma(k+1:end,1:k);
      
%% 3) Estimate the parameter zeta in equation (15) AER

zeta     = Sigmamu2*(Sigmamu1\eye(k)); % Sigmamu2*(Sigmamu1^{-1})

%% 4) Partition Sigma 

Sig11    = Sigma(1:k,1:k);      

Sig21    = Sigma(k+1:n,1:k); 

Sig22    = Sigma(k+1:n,k+1:n);

%% 5) Estimate the parameter eta in equation (16) AER

eta      = (Sig21' - (Sig11*zeta'))*((Sig22 - (Sig21*zeta'))\eye(n-k));

%% 6) Implement my version of equation (17) in Karel's AER

aux1     = (eye(k)-(eta*zeta))\eye(k);

aux      = [aux1 ; ...
           zeta*aux1];
       
%% 7) Estimate the parameter S1E[u_1t u_{1t'}] in equation 15)
% (under a Cholesky Assumption)

S1DS1    = Sig11 - (eta*Sig21) - (Sig21'*eta') + (eta*Sig22*eta');

if min(eig(Sigma))>0

S1Dchol  = chol(S1DS1, 'lower'); 

%% 8) Estimate the first column of B (up to scale)
            
B1aux       = aux*S1Dchol;

%% 9) Estimate the IRFs

B1        = x*B1aux(:,1)./B1aux(nvar,1);    
      
IRFSVARIV = reshape(sum(bsxfun(@times,Csim,B1'),2),[n,hori+1]);

else

IRFSVARIV = NaN(n,hori+1);    
    
end

end

