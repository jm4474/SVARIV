function  [bic,aic,hqic] = bicaic(series,pmax)
% -------------------------------------------------------------------------
% This function selects number of lags in the VAR model for "series"
% based on AIC,BIC and HQ information criteria
% 
% Inputs:
% - series: matrix of dimension (T times n) containing the time series
% - pmax: maximum possible number of lags  
% Outputs:
% - optimal lag-length by information criterion
% 
% This version: March 31, 2015
% Please, cite Gafarov, B. and Montiel-Olea, J.L. (2015) 
% "ON THE MAXIMUM AND MINIMUM RESPONSE TO AN IMPULSE IN SVARS"
% -------------------------------------------------------------------------
    

%% Definitions
[T,N]=size(series(pmax+1:end,:)); % consider first pmax periods as presample


%% Allocate memory
aic = zeros(pmax,1); 
bic = zeros(pmax,1); 
hqic = zeros(pmax,1);


%% Compute information criteria
for px = 1:pmax
    series_ = series(pmax-px+1:end,:);
    [~,~,eta,~] = RForm_VAR(series_,px); %Apply the function RForm_VAR.m to estimate reduced form parameters 
    lhd = log(det((eta*eta')/T));
    pty = px*N^2/T;
    aic(px)  = lhd + 2*pty;
    bic(px)  = lhd + log(T)*pty;
    hqic(px) = lhd + 2*log(log(T))*pty;
end


%% Optimal lag
[~,aic] = min(aic);
[~,bic] = min(bic); 
[~,hqic] = min(hqic); 

end
