function [reject, bootsIRFs] = GasydistbootsAR(ydata, T, seed, n, NB, p, norm, scale, horizons, confidence, SVARinp, NWlags, AL, Sigma, Gamma, V, WHatall, grid)
%  -Provides inference for SVAR-IV based on samples from the asy. dist.
%  -Syntax:
%    [IRFs, bootsIRFs] = Gasydistboots(seed, I, n, p, nvar, x, hori, confidence, vecAL, vechSigma, Gamma, Whatall)
%
%  -Inputs:
%     ydata: Endogenous variables from the VAR model
%         T: time periods
%      seed: seed structure  
%         n: number of variables in the VAR
%        NB: number of samples from the asymptotic distribution
%         p: number of lags in the VAR
%      norm: normalizing variable
%     scale: scale of the shock
%  horizons: number of horizons (IRFs) (does not include the impact or horizon 0)
%confidence: confidence level
%    NWlags: Newey-West lags
%        AL: point estimator of AL
%     Sigma: point estimator of Sigma
%     Gamma: point estimator of vec(Gamma)
%         V: Matrix such that vech(Sigma)=Vvec(Sigma)
%   Whatall: Covariance mat of (vec(AL)', vech(Sigma)', vec(Gamma)')
%      grid: Grid of values for the lambda to be used in the Anderson-Rubin Confidence Set
%
%  -Output:
%   reject: 4d logical array for whether an IRF is rejected or not, for
%   each lambda, variable, horizon, cumulative and non-cumulative
% bootsIRF: alpha/2 and 1-alpha/2 quantiles of IRFs
% 
% This version: July 17th, 2017
% Last edited by José Luis Montiel-Olea

%% 2) Set the seed and the main directory

rng(seed); clear seed

cd ..

cd ..

direct = pwd;

cd(direct)

%% 3) Create an RForm (if necessary)

if (nargin == 10)

        % This essentially checks whether the user provided or not his RForm. If
        % the user didn't, then we calculate it. If the user did, we skip this section and
        % use his/her RForm.

    SVARinp.ydata = ydata;

    SVARinp.Z = z;

    SVARinp.n        = size(ydata,2); %number of columns(variables)

    RForm.p          = p; %RForm_user.p is the number of lags in the model


    %a) Estimation of (AL, Sigma) and the reduced-form innovations


    % This essentially checks whether the user provided or not his RForm. If
    % the user didn't, then we calculate it. If the user did, we skip this section and
    % use his/her RForm.

    addpath('functions/RForm');

    [RForm.mu, ...
     RForm.AL, ...
     RForm.Sigma,...
     RForm.eta,...
     RForm.X,...
     RForm.Y]        = RForm_VAR(SVARinp.ydata, p);

    %RForm.AL(:),RForm.V*RForm.Sigma(:),RForm.Gamma(:), RForm.WHatall

    %b) Estimation of Gammahat (n times 1)

    RForm.Gamma      = RForm.eta*SVARinp.Z(p+1:end,1)/(size(RForm.eta,2));   %sum(u*z)/T. Used for the computation of impulse response.
    %(We need to take the instrument starting at period (p+1), because
    %we there are no reduced-form errors for the first p entries of Y.)

    %c) Add initial conditions and the external IV to the RForm structure

    RForm.Y0         = SVARinp.ydata(1:p,:);

    RForm.externalIV = SVARinp.Z(p+1:end,1);

    RForm.n          = SVARinp.n;

    %a) Covariance matrix for vec(A,Gammahat). Used
    %to conduct frequentist inference about the IRFs. 
    [RForm.WHatall,RForm.WHat,RForm.V] = ...
    CovAhat_Sigmahat_Gamma(p,RForm.X,SVARinp.Z(p+1:end,1),RForm.eta,NWlags);     

    %NOTES:
    %The matrix RForm.WHatall is the covariance matrix of 
    % vec(Ahat)',vech(Sigmahat)',Gamma')'

    %The matrix RForm.WHat is the covariance matrix of only
    % vec(Ahat)',Gamma')' 

    % The latter is all we need to conduct inference about the IRFs,
    % but the former is needed to conduct inference about FEVDs.
    
    AL = RForm.AL;
    
    Sigma = RForm.Sigma;
    
    Gamma = RForm.Gamma;
    
    V = RForm.V;
    
    WHatall = RForm.WHatall;
    
    n            = RForm.n;

elseif (nargin == 15)

    % We will use the user's RForm

else 

    % Error, not enough inputs.

end

vechSigma = V * Sigma(:);

AL = AL(:);

%% 4) Estimation of the asymptotic variance of A,Gamma

% Definitions


d            = ((n^2)*p)+(n);     %This is the size of (vec(A)',Gamma')'

%dall         = d+ (n*(n+1))/2;    %This is the size of (vec(A)',vec(Sigma), Gamma')'



%% 5) Make sure that Whatall is symmetric and positive semidefinite

%dall        = size(WHatall,1);  


dall        = size(WHatall,1);

WHatall     = (WHatall + WHatall')/2;
    
[aux1,aux2] = eig(WHatall);
    
WHatall     = aux1*max(aux2,0)*aux1'; 

%% 6) Generate draws
% Centered at (vec(AL)', Gamma')'


gvar    = [mvnrnd(zeros(NB,dall),(WHatall)/T)',...
                     zeros(dall,1)];
          %Added an extra column of zeros to access point estimators       
    
Draws   = bsxfun(@plus,gvar,...
          [AL;vechSigma;Gamma(:)]);

k       = size(Gamma,1)/n;        

%% 7 Evaluate the parameter of interest  

ndraws        = size(Draws,2);
     
pdSigma       = zeros(1,ndraws);

addpath('functions/StructuralIRF');

RFormIRFBoots = zeros(n, horizons + 1,ndraws,2); %4th dimension corresponds to non-cumulative and cumulative values.

AlphaBoots    = zeros(1, ndraws);

for idraws    = 1:ndraws
    
    %i) Generate the draws for AL 
        
    AL        = reshape(Draws(1:(n^2)*p,idraws),[n,n*p]);

    %ii) Draws from Sigma

    vechSigma = Draws((n^2)*p+1:(n^2)*p+(n*(n+1)/2),idraws);

    Sigma     = tril(ones(n),0); Sigma(Sigma==1) = vechSigma';

    Sigma     = Sigma + tril(Sigma,-1)';

    %Check if the draws are positive definite

    if min(eig(Sigma))>0

        pdSigma(1,idraws) = 1;

    else

        pdSigma(1,idraws) = 0;

    end

    %iii) Draws from Gamma

    Gamma = reshape(Draws(((n^2)*p)+(n*(n+1)/2)+1:end,idraws),...
      [n,k]);  

    [RFormIRFBoots(:,:,idraws,:), AlphaBoots(idraws)] = IRFSVARIV(AL,Sigma,Gamma,horizons,scale,norm);
    
    RFormIRFBoots(:,:,idraws,:) = RFormIRFBoots(:,:,idraws,:).*AlphaBoots(idraws);
    % RFormIRFBoots(n,horizons+1,idraws,2)
    
    clear AL vechSigma Gamma 
    
end

grid_size       = size(grid,3);

test_aux      = zeros(grid_size, ndraws, n, horizons+1,2); % 5th dimension is for cumulative and non-cumulative

for var         = 1:n

    for horizon = 1:horizons+1
        
        aux_grid = reshape(grid(var,horizon,:),[grid_size,1]);
        
        test_aux(:,:,var,horizon,:) = TestStatistic(var, horizon, RFormIRFBoots, AlphaBoots, aux_grid, T);
        
    end
    
end

AR_test         = (test_aux - test_aux(:,ndraws,:,:,:)).^2;
% AR_test(grid_size, ndraws, n, horizons+1,cumulative)

%% 8) Implement "Standard" Bootstrap Inference

aux          = reshape(pdSigma,[1,ndraws]);


bootsIRFs    = quantile(AR_test(:,aux==1,:,:,:),...
                          [((1-confidence)/2),1-((1-confidence)/2)],2);      
% bootsIRFs(grid_size, confidence intervals, variables, horizons+1, cumulative)

                      
test_T = test_aux(:,1001,:,:,:);
% differnece_T(grid_size, ndraws, n, horizons+1,cumulative)

% check whether each one would be rejected or not

reject       = (test_T(:,1,:,:,:) < bootsIRFs(:,1,:,:,:)) | (test_T(:,1,:,:,:) > bootsIRFs(:,2,:,:,:));

reject       = reshape(reject,[grid_size, n, horizons+1,2]);
% reject(grid_size, variables, horizons+1, cumulative)








