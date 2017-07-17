function [C] = MARep(AL,p,hori)
% -------------------------------------------------------------------------
% Transforms the A(L) parameters of a reduced-form VAR
% into the coefficients C of the MA representation.
% 
% Inputs:
% - AL: VAR model coefficients
% - p: number of lags in the VAR model
% - hori: forecast horizon
% Outputs:
% - C: MA representation coefficients
% 
% This version: July 17, 2017
% -------------------------------------------------------------------------

%% 1) Reshape AL into a 3-D array

n         = size(AL,1);

vecAL     = reshape(AL,[n,n,p]); 

%% 2) Initialize the value of the auxiliary array vecALrevT

vecALrevT = zeros(n,n,hori);

for ihori     = 1:hori
    
    if  ihori < (hori-p)+1
        
        vecALrevT(:,:,ihori) = zeros(n,n);
    
    else
        
        vecALrevT(:,:,ihori) = vecAL(:,:,(hori-ihori)+1)';
    end
    
end

vecALrevT     = reshape(vecALrevT,[n,n*hori]);

% vecALrevT   = [0,0, ... vecAL(:,:,p), ... vecAL(:,:,1)];

%% MA coefficients

C             = repmat(vecAL(:,:,1),[1,hori]);

for     ihori = 1:hori-1
    
    C(:,(n*ihori)+1:(n*(ihori+1))) = [eye(n),C(:,1:n*ihori)] ...
        * vecALrevT(:,(hori*n-(n*(ihori+1)))+1:end)';
    
end

end

