function GasydistbootsAR_check(ydata, T, n, NB, p, norm, scale, horizons, confidence, NWlags, AL, Sigma, Gamma, V, WHatall, multiplier, grid_size, ARylim)

% Checks whether the inputs from Bootstrap_Plots.m are valid.
%  -Inputs:
%     ydata: Endogenous variables from the VAR model
%         T: time periods
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
%multiplier: Scalar that GasydistbootsAR will use to create the "grid" of null hypothesis for IRFs
% grid_size: Size of the grid of values for the lambda to be used in the Anderson-Rubin Confidence Set

% This version: Semptember 14th, 2018
% Last edited by José Luis Montiel-Olea

%% Check inputs

%check ydata
if isempty(ydata) == 1
    
    error('ydata must be a (T times n) matrix. It is now empty');

end

if isnumeric(ydata) == 0
    
    error('ydata must be numeric.');
    
end

if ismatrix(ydata) == 0
    
    error('ydata must be two dimensional (T times n).');

end 

if size(ydata,1) < size(ydata,2)
    
    warning('ydata: number of rows < number of columns. ydata must be T times n.');
    
end

if size(ydata,1) ~= T
    
    error('ydata must have T rows');
    
end

%check T
if isempty(T)
    
    error('T must be assigned a value.');

end

if isnumeric(T) == 0
    
    error('T must be numeric.');

end

if numel(T) ~= 1
    
    error('T must have only one element');
    
end

%Check n
if isempty(n)
    
    error('n must be assigned a value.');

end

if isnumeric(n) == 0
    
    error('n must be numeric.');

end

if numel(n) ~= 1
    
    error('n must have only one element');
    
end

if floor(n) ~= n
    
    error('n must be an integer.');
    
end

%Check NB
if isempty(NB)
    
    error('NB must be assigned a value.');

end

if isnumeric(NB) == 0
    
    error('NB must be numeric.');

end

if numel(NB) ~= 1
    
    error('NB must have only one element');
    
end

if floor(NB) ~= NB
    
    error('NB must be an integer.');
    
end

%Check p
if isempty(p)
    
    error('p must be assigned a value.');

end

if isnumeric(p) == 0
    
    error('p must be numeric.');

end

if numel(p) ~= 1
    
    error('p must have only one element');
    
end

if floor(p) ~= p
    
    error('p must be an integer.');
    
end

%check norm
if isempty(norm)
    
    error('norm must be assigned a value.');

end

if isnumeric(norm) == 0
    
    error('norm must be numeric.');

end

if numel(norm) ~= 1
    
    error('norm must have only one element');
    
end

if floor(norm) ~= norm
    
    error('norm must be an integer.');
    
end

if norm > n
    
    error('norm must be smaller or equal to the total amount of variables');
    
end

if norm < 1
    
    error('norm must be >= 1');
    
end

%check scale
if isempty(scale)
    
    error('scale must be assigned a value.');

end

if isnumeric(scale) == 0
    
    error('scale must be numeric.');

end

if numel(scale) ~= 1
    
    error('scale must have only one element');
    
end

%check horizons
if isempty(horizons)
    
    error('horizons must be assigned a value.');

end

if isnumeric(horizons) == 0
    
    error('horizons must be numeric.');

end

if numel(horizons) ~= 1
    
    error('horizons must have only one element');
    
end

if horizons < 0
    
    error('horizons must be >= 0.');
    
end

if floor(horizons) ~= horizons
    
    error('horizons must be an integer.');
    
end

%check confidence
if isempty(confidence)
    
    error('confidence must be assigned a value.');

end

if isnumeric(confidence) == 0
    
    error('confidence must be numeric.');

end

if numel(confidence) ~= 1
    
    error('confidence must have only one element');
    
end

if(confidence <= 0 || confidence >= 1)
    
    error('confidence must be > 0 and < 1');
    
end

%check NWlags
if isempty(NWlags)
    
    error('NWlags must be assigned a value.');

end

if isnumeric(NWlags) == 0
    
    error('NWlags must be numeric.');

end

if numel(NWlags) ~= 1
    
    error('NWlags must have only one element');
    
end

if NWlags < 0
    
    error('NWlags must be positive');
    
end

%check AL
if isempty(AL) == 1
    
    error('AL must be a (n times (n*p)) matrix. It is now empty');

end

if isnumeric(AL) == 0
    
    error('AL must be numeric.');
    
end

if ismatrix(AL) == 0
    
    error('AL must be two dimensional (n times (n*p)).');

end 

if size(AL,1) ~= n
    
    error('size(AL,1) must be n') ;
    
end

if size(AL,2) ~= n*p
    
    error('size(AL,2) must be n*p');
    
end

%check Sigma
if isempty(Sigma) == 1
    
    error('Sigma must be a (n times n) matrix. It is now empty');

end

if isnumeric(Sigma) == 0
    
    error('Sigma must be numeric.');
    
end

if ismatrix(Sigma) == 0
    
    error('Sigma must be two dimensional (n times n).');

end 

if size(Sigma,1) ~= n
    
    error('size(Sigma,1) must be n') ;
    
end

if size(Sigma,2) ~= n
    
    error('size(Sigma,2) must be n');
    
end

%check Gamma
if isempty(Gamma) == 1
    
    error('Gamma must be a (n times n) matrix. It is now empty');

end

if isnumeric(Gamma) == 0
    
    error('Gamma must be numeric.');
    
end

if isvector(Gamma) == 0
    
    error('Gamma must be an array (n times 1).');

end 

if size(Gamma,1) ~= n
    
    error('size(Gamma,1) must be n') ;
    
end

%check V
if isempty(V) == 1
    
    error('V must be a (factorial(n) times n*n) matrix. It is now empty');

end

if isnumeric(V) == 0
    
    error('V must be numeric.');
    
end

if ismatrix(V) == 0
    
    error('V must be a matrix.');

end 

if size(V,1) ~= factorial(n)
    
    error('size(V,1) must be factorial(n)') ;
    
end

if size(V,2) ~= n*n
    
    error('size(V,2) must be n*n');
    
end

%check WHatall
%WHatall is the covariance matrix of vec(A),vech(Sigma),Gamma
if isempty(WHatall) == 1
    
    error('WHatall must be a matrix. It is now empty');

end

if isnumeric(WHatall) == 0
    
    error('WHatall must be numeric.');
    
end

if ismatrix(WHatall) == 0
    
    error('WHatall must be a matrix.');

end 

%check multiplier

if isempty(multiplier) == 1
    
    error('multiplier must be a scalar (1 times 1). It is now empty');

end

if isnumeric(multiplier) == 0
    
    error('multiplier must be numeric.');
    
end

if numel(multiplier) ~= 1
    
    error('multiplier must have only one element');
    
end

if multiplier <= 0
    
    error('multiplier must be > 0');
    
end

%check grid_size

if isempty(grid_size) == 1
    
    error('grid_size must be a scalar (1 times 1). It is now empty');

end

if isnumeric(grid_size) == 0
    
    error('grid_size must be numeric.');
    
end

if numel(grid_size) ~= 1
    
    error('grid_size must have only one element');
    
end

if grid_size <= 0
    
    error('grid_size must be > 0');
    
end

%check ARylim

if isempty(ARylim) == 1
    
    error('ARylim must be a three dimensional array (n times 2 times 2). It is now empty');

end

if isnumeric(ARylim) == 0
    
    error('ARylim must be numeric.');
    
end

if size(ARylim,1) ~= n
    
    error('ARylim''s first row size must be n');
    
end

if size(ARylim,2) ~= 2
    
    error('ARylim''s second row size must be 2');
    
end

if size(ARylim,3) ~= 2
    
    error('ARylim''s third row size must be 2');
    
end

end
    