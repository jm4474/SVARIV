function Gasydistboots_Check(NB, n, p, norm, scale, horizons, confidence, T, f, NWlags, AL, Sigma, Gamma, V, WHatall)
% Checks whether the inputs from SVARIV.m are valid.
%  -Syntax:
%    Gasydistboots(seed, I, n, p, nvar, x, hori, confidence, vecAL, vechSigma, Gamma, Whatall)
%  -Inputs:
%       NB: number of samples from the asymptotic distribution
%        n: number of variables in the VAR
%        p: number of lags in the VAR
%     norm: normalizing variable
%    scale: scale of the shock
% horizons: number of horizons (IRFs) (does not include the impact or horizon 0)
%confidence: confidence level
%        T: time periods
%        f: function handle (depends on AL, Sigma, Gamma, hori, x, nvar)
%   NWlags: Newey-West lags
%       AL: point estimator of AL
%vechSigma: point estimator of vech(Sigma)
%    Gamma: point estimator of vec(Gamma)
%  Whatall: Covariance mat of (vec(AL)', vech(Sigma)', vec(Gamma)')     

% This version: Semptember 14th, 2018
% Last edited by José Luis Montiel-Olea

%% 1) Check if the inputs are valid

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

%check f
if isempty(f) == 1
    
    error('f must not be empty');
    
end

if isa(f, 'function_handle') == 0
    
    error('f must be a function handle');
    
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

end