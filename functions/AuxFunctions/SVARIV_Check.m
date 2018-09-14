function SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name)
% Checks whether the inputs from SVARIV.m are valid.
%-Syntax:
%       SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name)
%
% -Inputs:
%       p:            Number of lags in the VAR model                                                    (1 times 1)                                          
%       confidence:   Value for the standard and weak-IV robust confidence set                           (1 times 1) 
%       ydata:        Endogenous variables from the VAR model                                            (T times n) 
%       z:            External instrumental variable                                                     (T times 1)
%       NWlags:       Newey-West lags                                                                    (1 times 1)
%       norm:         Variable used for normalization                                                    (1 times 1)
%       scale:        Scale of the shock                                                                 (1 times 1)
%       horizons:     Number of horizons for the Impulse Response Functions (IRFs) 
%                     (does not include the impact horizon 0)                                            (1 times 1)
%       savdir:       Directory where the figures generated will be saved                                (String)
%       columnnames:  Vector with the names for the endogenous variables, in the same order as ydata     (1 times n)
%       IRFselect:    Indices for the variables that the user wants separate IRF plots for               (1 times q)
%       cumselect:    Indices for the variables that the user wants cumulative IRF plots for
%       time:         Time unit for the dataset (e.g. year, month, etc.)                                 (String)
%       dataset_name: The name of the dataset used for generating the figures (used in the output label) (String)

% This version: Semptember 14th, 2018
% Last edited by José Luis Montiel-Olea
%% 1) Check if the inputs are valid

if ~strncmp(version, '9', 1)
    
    warning('You are using a version other than v9, in which this program was tested.')

end 

%check p value
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

%check confidence value
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

T = size(ydata,1);

%check z
if isempty(z) == 1
    
    error('z must be (T times 1) vector. It is now empty');

end

if isnumeric(z) == 0
    
    error('z must be numeric');
    
end

if ismatrix(z) == 0
    
    error('z must be two dimensional (T times 1)');
    
end

if size(z,2) ~= 1
    
    error('There must be only one instrument (z must have only one column)');
    
end

if size(z,1) ~= T
    
    error('z must have T rows');
    
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

if floor(NWlags) ~= NWlags
    
    error('NWlags must be an integer.');
    
end

%Check norm
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

if norm > size(ydata, 2)
    
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

%check savdir
if isempty(savdir)
    
    savdir('savdir must be assigned a value.');

end

if ischar(savdir) == 0
    
    error('savdir must be a character array.');

end

if size(savdir,1) > 1
    
    error('savdir must only contain one string');
    
end

%columnnames
if isempty(columnnames)
    
    error('columnnames must be assigned a character array.');

end

if iscellstr(columnnames) == 0
    
    error('columnnames must be a cell array of character vectors.');

end

if size(columnnames,2) > size(ydata, 2)
    
     error('columnames size must be smaller or equal to the total amount of variables');
     
end

%Check IRFselect
if isempty(IRFselect) ~= 0
   
    %IRFselect can be empty. If it isn't, then we need to check if it
    %follows the following criteria:
    if isnumeric(IRFselect) == 0
        
        error('IRFselect must be a numeric array');
        
    end
    
    if ismatrix(IRFselect) ~= 2
        
        error('IRFselect must be a row vector');
    
    end
    
    if size(IRFselect,1) > 1
        
        error('IRFselect should be a row vector');
        
    end
    
    if size(IRFselect,2) > size(ydata, 2)
        
        error('IRFselect size must be smaller or equal to the total amount of variables');
        
    end
    
    
    
end

%Check cumselect
if isempty(cumselect) ~= 0
   
    %cumselect can be empty. If it isn't, then we need to check if it
    %follows the following criteria:
    if isnumeric(cumselect) == 0
        
        error('cumselect must be a numeric array');
        
    end
    
    if ismatrix(cumselect) ~= 2
        
        error('cumselect must be a row vector');
    
    end
    
    if size(cumselect,1) > 1
        
        error('cumselect must be a row vector');
        
    end
    
    if size(cumselect,2) > size(ydata, 2)
        
        error('cumselect size must be smaller or equal to the total amount of variables');
        
    end
    
end
    
%check time

if ischar(time) == 0
    
    error('time must be a character array');
    
end

%we need to check whether time's size is 1x(number of characters).
if ismatrix(time) ~= 1
    
    error('time must be a character array.');
    
end

if size(time,1) ~= 1
    
    error('time must be a character array.');
    
end

%check dataset_name
    
if ischar(dataset_name) == 0
    
    error('dataset_name must be a character array');
    
end

%we need to check whether time's size is 1x(number of characters).
if ismatrix(time) ~= 1
    
    error('dataset_name must be a character array.');
    
end

if size(time,1) ~= 1
    
    error('dataset_name must be a character array.');
    
end

end

    



    


 
    
    





