function [InferenceMSW,Plugin,Chol] = MSWfunction(confidence,nvar,scale,horizons,RForm,display_on)
% -Reports the confidence interval for IRF coefficients described in Montiel-Olea, Stock, and Watson (2017). This version: July 11th, 2017. 
% -Syntax:
%       [InferenceMSW,Plugin,Chol] = MSWfunction(confidence,nvar,x,hori,RForm,display_on)
% -Inputs:
%    confidence: confidence level                            (1 x 1)
%          nvar: variable defining the normalization         (1 x 1)
%         scale: scale of the shock                          (1 x 1)
%      horizons: Number of horizons for the Impulse Response 
%                Functions (IRFs) (does not include the 
%                impact or horizon 0)                        (1 x 1)
%         RForm: reduced-form structure.                     (structure)
%                
%                The structure must contain the following fields:
%                AL= n x np matrix of reduced-form coefficients
%                p = number of lags used in the estimation of AL
%                n = dimension of the SVAR
%            Sigma = n x n matrix of reduced-form residuals
%       RForm.What = (n^2p + n) covariance matrix of vec(A),Gamma
%        RForm.eta = n x T matrix of reduced-form residuals
%      RForm.Gamma = Estimator of E[z_t eta_t]

%    display_on: dummy variable. 
%
% -Output:
%  InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%        Plugin: Structure containing standard plug-in inference
%          Chol: Cholesky IRFs
%

critval = norminv(1-((1-confidence)/2),0,1)^2;

%% 1) Create the MA coefficients based on the AL matrix

addpath('functions/StructuralIRF') 

Caux        = [eye(RForm.n),MARep(RForm.AL,RForm.p,horizons)]; 
%The function MARep uses the reduced-form estimator vec(A) to 
%compute the implied MA coefficients. You can replace this function
%by your own routine, but make sure that the dimensions match.
%Using the data from SVARIV_tax, Caux is a 9x63 matrix. 

C           = reshape(Caux,[RForm.n,RForm.n,horizons+1]); 
%Put Caux in a 3-dimensional array
%Each page in the third dimension represents the MA coefficients for that
%time horizon

Ccum        = cumsum(C,3); 
%Compute the cumulative MA coefficients as we report both 
%cumulative and noncumulative IRFs

[G,Gcum]    = Gmatrices(RForm.AL,MARep(RForm.AL,RForm.p,horizons),RForm.p,horizons,RForm.n); 
%(and the derivatives).

%% 2) Compute the Cholesky estimator for comparison

B1chol      = chol(RForm.Sigma)'; 
%Compute the Cholesky estimator for comparison

B1chol      = scale*(B1chol(:,1)./B1chol(nvar,1));

Chol(:,:,1) = reshape(sum(bsxfun(@times,C,B1chol'),2),[RForm.n,horizons+1]);

Chol(:,:,2) = reshape(sum(bsxfun(@times,Ccum,B1chol'),2),[RForm.n,horizons+1]);

%% 3) Label the submatrices of the asy var of (vecA,Gamma)  

W1          = RForm.WHat(1:(RForm.n^2)*RForm.p,1:(RForm.n^2)*RForm.p);

W12         = RForm.WHat(1:(RForm.n^2)*RForm.p,1+(RForm.n^2)*RForm.p:end);

W2          = RForm.WHat(1+(RForm.n^2)*RForm.p:end,1+(RForm.n^2)*RForm.p:end);

%Up to here, most of the computations were standard. From now on, we 
%will rely heavily on the notation in Montiel-Stock-Watson. 

%% 4) Definitions to apply the formulae in MSW for noncumulative IRFs

%a) Definitions to compute the MSW confidence interval for $\lambda_{k,i}$

n         = RForm.n;  %dimension of the SVAR

T         = (size(RForm.eta,2));  % Number of observations (number of time periods)

e         = eye(n); 

ahat      = zeros(n,horizons+1); 

bhat      = zeros(n,horizons+1);

chat      = zeros(n,horizons+1);

Deltahat  = zeros(n,horizons+1);
    
MSWlbound = zeros(n,horizons+1);

MSWubound = zeros(n,horizons+1);

casedummy = zeros(n,horizons+1);

for j =1:n
    
    for ih=1:horizons+1
        
    ahat(j,ih)     = (T*(RForm.Gamma(nvar,1)^2))-(critval*W2(nvar,nvar));
    
    bhat(j,ih)     = -2*T*scale*(e(:,j)'*C(:,:,ih)*RForm.Gamma)*RForm.Gamma(nvar,1)...
        + 2*critval*scale*(kron(RForm.Gamma',e(:,j)'))*G(:,:,ih)*W12(:,nvar)...
        + 2*critval*scale*e(:,j)'*C(:,:,ih)*W2(:,nvar);
    
    chat(j,ih)     = ((T^.5)*scale*e(:,j)'*C(:,:,ih)*RForm.Gamma).^2 ...
        -critval*(scale^2)*(kron(RForm.Gamma',e(:,j)'))*G(:,:,ih)*W1*...
        ((kron(RForm.Gamma',e(:,j)'))*G(:,:,ih))' ...
        -2*critval*(scale^2)*(kron(RForm.Gamma',e(:,j)'))*G(:,:,ih)*W12*C(:,:,ih)'*e(:,j)...
        -critval*(scale^2)*e(:,j)'*C(:,:,ih)*W2*C(:,:,ih)'*e(:,j); 
    
    Deltahat(j,ih) = bhat(j,ih).^2-(4*ahat(j,ih)*chat(j,ih));
    
    if ahat(j,ih)>0 && Deltahat(j,ih)>0
        
   casedummy(j,ih) = 1;
   
   MSWlbound(j,ih) = (-bhat(j,ih) - (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
        
   MSWubound(j,ih) = (-bhat(j,ih) + (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
        
    elseif ahat(j,ih)<0 && Deltahat(j,ih)>0
        
   casedummy(j,ih) = 2;
   
   MSWlbound(j,ih) = (-bhat(j,ih) + (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
   
   MSWubound(j,ih) = (-bhat(j,ih) - (Deltahat(j,ih)^.5))/(2*ahat(j,ih));
   
    elseif ahat(j,ih)>0 && Deltahat(j,ih)<0
        
   casedummy(j,ih) = 3;
   
   MSWlbound(j,ih) = NaN;
   
   MSWubound(j,ih) = NaN;
   
    else 
        
   casedummy(j,ih) = 4;
   
   MSWlbound(j,ih) = -inf;
   
   MSWubound(j,ih) = inf;
   
    end
    
    end
end

MSWlbound(nvar,1)=scale;
    
MSWubound(nvar,1)=scale;
    
%% 5) Save all the output in the structure Inference.MSW    
    
InferenceMSW.ahat      = ahat; clear ahat
    
InferenceMSW.bhat      = bhat; clear bhat

InferenceMSW.chat      = chat; clear chat

InferenceMSW.Deltahat  = Deltahat; clear Deltahat

InferenceMSW.casedummy = casedummy; clear casedummy

InferenceMSW.MSWlbound = MSWlbound; clear MSWlbound

InferenceMSW.MSWubound = MSWubound; clear MSWubound

InferenceMSW.T         = T;

%MSWL bound and MSWubound contain the bounds of the confidence interval
%These bounds can be infinity depending on the cases described in
%the paper. Note that in case 2, we make a slight abuse of notation
%as MSWlbound in fact refers to the largest value in the first
%open ray (-infty, c) that defines the MSW confidence interval. 

%% 6) For comparison purposes, we report the standard delta-method CI    

lambdahat                  = zeros(n,horizons+1);

DmethodVar                 = zeros(n,horizons+1);

Dmethodlbound              = zeros(n,horizons+1);

Dmethodubound              = zeros(n,horizons+1);

for ih = 1:horizons + 1
    
    for ivar = 1:n
        
        lambdahat(ivar,ih)     = scale*e(:,ivar)'*C(:,:,ih)*RForm.Gamma./RForm.Gamma(nvar,1);
        
        d1                     = (kron(RForm.Gamma',e(:,ivar)')*scale*G(:,:,ih));
        
        d2                     = (scale*e(:,ivar)'*C(:,:,ih))-(lambdahat(ivar,ih)*e(:,nvar)');                                    
         
        d                      = [d1,d2]';         
        
        DmethodVar(ivar,ih)    = d'*RForm.WHat*d;           
        
        Dmethodlbound(ivar,ih) = lambdahat(ivar,ih)-...
            ((critval./T)^.5)*(DmethodVar(ivar,ih)^.5)/abs(RForm.Gamma(nvar,1));
        
        Dmethodubound(ivar,ih) = lambdahat(ivar,ih)+...
            ((critval./T)^.5)*(DmethodVar(ivar,ih)^.5)/abs(RForm.Gamma(nvar,1));  
        
        clear d1 d2 d;
    
    end
end

%% 7) Save the delta-method output in Inference.MSW and also in 
% the Plug-in structure

InferenceMSW.Dmethodlbound = Dmethodlbound;                                        clear Dmethodlbound;

InferenceMSW.Dmethodubound = Dmethodubound;                                        clear Dmethodubound;    

Plugin.IRF                 = lambdahat;                                            clear lambdahat; 

Plugin.IRFstderror         = (DmethodVar.^.5)./((T^.5)*abs(RForm.Gamma(nvar,1)));  clear DmethodVar;

%% 8) Definitions to apply the formulae in MSW for cumulative IRFs     

ahatcum      = zeros(n,horizons+1);

bhatcum      = zeros(n,horizons+1);

chatcum      = zeros(n,horizons+1);

Deltahatcum  = zeros(n,horizons+1);

MSWlboundcum = zeros(n,horizons+1);

MSWuboundcum = zeros(n,horizons+1);

casedummycum = zeros(n,horizons+1);

for j = 1:n
    for ih = 1:horizons+1

        ahatcum(j,ih)          = (T*(RForm.Gamma(nvar,1)^2))-(critval*W2(nvar,nvar));

        bhatcum(j,ih)          = -2*T*scale*(e(:,j)'*Ccum(:,:,ih)*RForm.Gamma)*RForm.Gamma(nvar,1)...
            + 2*critval*scale*(kron(RForm.Gamma',e(:,j)'))*Gcum(:,:,ih)*W12(:,nvar)...
            + 2*critval*scale*e(:,j)'*Ccum(:,:,ih)*W2(:,nvar);

        chatcum(j,ih)          = ((T^.5)*scale*e(:,j)'*Ccum(:,:,ih)*RForm.Gamma).^2 ...
            -critval*(scale^2)*(kron(RForm.Gamma',e(:,j)'))*Gcum(:,:,ih)*W1*...
            ((kron(RForm.Gamma',e(:,j)'))*Gcum(:,:,ih))' ...
            -2*critval*(scale^2)*(kron(RForm.Gamma',e(:,j)'))*Gcum(:,:,ih)*W12*Ccum(:,:,ih)'*e(:,j)...
            -critval*(scale^2)*e(:,j)'*Ccum(:,:,ih)*W2*Ccum(:,:,ih)'*e(:,j); 

        Deltahatcum(j,ih)      = (bhatcum(j,ih).^2)-(4*ahatcum(j,ih)*chatcum(j,ih));

        if (ahatcum(j,ih)>0 && Deltahatcum(j,ih)>0)

            casedummycum(j,ih) = 1;

            MSWlboundcum(j,ih) = (-bhatcum(j,ih) - (Deltahatcum(j,ih)^.5))/(2*ahatcum(j,ih));

            MSWuboundcum(j,ih) = (-bhatcum(j,ih) + (Deltahatcum(j,ih)^.5))/(2*ahatcum(j,ih));

        elseif (ahatcum(j,ih)<0 && Deltahatcum(j,ih)>0)

            casedummycum(j,ih) = 2;

            MSWlboundcum(j,ih) = (-bhatcum(j,ih) + (Deltahatcum(j,ih)^.5))/(2*ahatcum(j,ih));

            MSWuboundcum(j,ih) = (-bhatcum(j,ih) - (Deltahatcum(j,ih)^.5))/(2*ahatcum(j,ih));

        elseif (ahatcum(j,ih)>0 && Deltahatcum(j,ih)<0)

            casedummycum(j,ih) = 3;

            MSWlboundcum(j,ih) = NaN;

            MSWuboundcum(j,ih) = NaN;

        else

            casedummycum(j,ih) = 4;

            MSWlboundcum(j,ih) = -inf;

            MSWuboundcum(j,ih) = inf;

        end
    end
end

MSWlboundcum(nvar,1) = scale;

MSWuboundcum(nvar,1) = scale;

%% 9) Save all the output in the structure Inference.MSW    

InferenceMSW.ahatcum      = ahatcum;      clear ahatcum

InferenceMSW.bhatcum      = bhatcum;      clear bhatcum

InferenceMSW.chatcum      = chatcum;      clear chatcum

InferenceMSW.Deltahatcum  = Deltahatcum;  clear Deltahatcum

InferenceMSW.MSWlboundcum = MSWlboundcum; clear MSWlboundcum

InferenceMSW.MSWuboundcum = MSWuboundcum; clear MSWuboundcum

InferenceMSW.casedummycum = casedummycum; clear casedummy
    
%% 10) We also report the standard delta-method CI for cumulative IRFs   

lambdahatcum     = zeros(n,horizons+1);

DmethodVarcum    = zeros(n,horizons+1);

Dmethodlboundcum = zeros(n,horizons+1);

Dmethoduboundcum = zeros(n,horizons+1);
        
for ih = 1:horizons + 1
    
    for ivar = 1:n
        
        lambdahatcum(ivar,ih)     = scale*e(:,ivar)'*Ccum(:,:,ih)*RForm.Gamma./RForm.Gamma(nvar,1);
        
        d1                        = (kron(RForm.Gamma',e(:,ivar)')*scale*Gcum(:,:,ih));
        
        d2                        = scale*e(:,ivar)'*Ccum(:,:,ih)-lambdahatcum(ivar,ih)*e(:,nvar)';                                    
        
        d                         = [d1,d2]';         
        
        DmethodVarcum(ivar,ih)    = d'*RForm.WHat*d;           
        
        Dmethodlboundcum(ivar,ih) = lambdahatcum(ivar,ih)-...
            ((critval./T)^.5)*(DmethodVarcum(ivar,ih)^.5)/abs(RForm.Gamma(nvar,1));
        
        Dmethoduboundcum(ivar,ih) = lambdahatcum(ivar,ih)+...
            ((critval./T)^.5)*(DmethodVarcum(ivar,ih)^.5)/abs(RForm.Gamma(nvar,1));  
        
        clear d1 d2 d;
    
    end
    
end

%% 11) Save the delta-method output in Inference.MSW and also in 
% the Plug-in structure

InferenceMSW.Dmethodlboundcum = Dmethodlboundcum;                                        clear Dmethodlboundcum;

InferenceMSW.Dmethoduboundcum = Dmethoduboundcum;                                        clear Dmethoduboundcum;

Plugin.IRFcum                 = lambdahatcum;                                            clear lambdahatcum; 

Plugin.IRFstderrorcum         = (DmethodVarcum.^.5)./((T^.5)*abs(RForm.Gamma(nvar,1)));  clear DmethodVarcum;
    
%% 12) Display the information about the estimated confidence interval
% and the first-stage. This part of the program is Optional.
if display_on==1
    
    disp('-');
    
    disp('NOTE: The Wald statistic for the covariance between the instrument and the normalized variable is:')
    
    Waldstat              = (((T^.5)*RForm.Gamma(nvar,1))^2)/RForm.WHat(((n^2)*RForm.p)+nvar,((n^2)*RForm.p)+nvar);
    
    display(Waldstat);
    
    InferenceMSW.Waldstat = Waldstat;
    
    fprintf('Given the confidence level, if the Wald statistic is larger than: \n\n');
    
    disp(critval);
    
    disp('The weak-IV robust confidence set will be a bounded interval for every horizon (check "casedummy" if not).')
    
    disp('--');
    
    %display('Also, the Wald statistic for the covariance between the instrument and the full vector of reduced-form residuals is')
    %WaldstatFull= (T)*(RForm.Gamma'*(RForm.WHat(((n^2)*RForm.p)+nvar:end,((n^2)*RForm.p)+nvar:end))^(-1)*RForm.Gamma);
    %display(WaldstatFull);
    %display('The 1-alpha quantile of this statistic is:')
    %display(chi2inv(confidence,n));
else
    
end

%% 13) Finally, we report the plug-in estimators of the target shock,
% the correlation between the structural shock and the instrument, and
% the forecast error variance decomposition. 

%9.1) Target Structural Shock:

%e) Estimated shock:
    
    Plugin.epsilonhat    = scale*RForm.Gamma'*(RForm.Sigma^(-1))*RForm.eta./RForm.Gamma(nvar,1);
    
    Plugin.epsilonhatstd = (Plugin.epsilonhat-mean(Plugin.epsilonhat))./std(Plugin.epsilonhat);
    
end

