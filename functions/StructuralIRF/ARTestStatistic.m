function [test] = ARTestStatistic(var, horizon, RFormIRFBoots, AlphaBoots, null_vec, T, ndraws)
%  -Provides inference for SVAR-IV based on samples from the asy. dist.
%  -Syntax:
%    [test] = TestStatistic(var, horizon, RFormIRFBoots, AlphaBoots, grid, T)
%  -Inputs:
%     var           : Variable for which we are running the test statistic
%     horizon       : Horizon for which we are running the test statistic
%     RFormIRFBoots : RForm impulse reponses
%     AlphaBoots    : Gama(1,1)
%     grid          : Grid of lambdas for each variable and horizon.
%     T             : Number of time periods
%
%  -Output:
%     test          : test statistic used in the bootstrap implementation
%                     of the Anderson-Rubin Confidence Set
% 
% This version: Semptember 4th, 2018
% Last edited by José Luis Montiel-Olea

%% 1) Main definitions and function

IRFBootsVH = RFormIRFBoots(var,horizon,:,:);
        
IRFBootsVH = reshape(IRFBootsVH, [1, ndraws,2]);

test(:,:,1) = (IRFBootsVH(:,:,1) - (null_vec * AlphaBoots(1,:)))*(T^.5); 
% (scale*Ck(A)*Gamma - lambda*Gamma(1,1))*T^(1/2)
% The third dimension allows for cumulative and non-cumulative computation.
        
test(:,:,2) = (IRFBootsVH(:,:,2) - (null_vec * AlphaBoots(1,:)))*(T^.5);

aux = test(:,ndraws,:);

test = test(:,1:(ndraws-1),:) - test(:,ndraws,:);

test(:,ndraws,:) = aux;
        
end