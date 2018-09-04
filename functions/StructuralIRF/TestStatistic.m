function [test] = TestStatistic(var, horizon, RFormIRFBoots, AlphaBoots, grid, T)
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
        
IRFBootsVH = reshape(IRFBootsVH, 1, 1001,2);

test(:,:,1) = (IRFBootsVH(:,:,1) - (grid * AlphaBoots(1,:)))*(T^.5); 
% (scale*Ck(A) - lambda)*Gamma*T^(1/2)
% The third dimension allows for cumulative and non-cumulative computation.
        
test(:,:,2) = (IRFBootsVH(:,:,2) - (grid * AlphaBoots(1,:)))*(T^.5);
        
end