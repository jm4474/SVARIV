function Bootstrap_Plots_Check(n,p,horizons,confidence, figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect)
% Checks whether the inputs from Bootstrap_Plots.m are valid.
%  -Syntax:
%    Bootstrap_Plots_Check(n,p,horizons,confidence, figureorder,Plugin,InferenceMSW,time,columnnames,savdir,direct,dataset_name,IRFselect,cumselect)
%  -Inputs:
%      n:               Number of variables in the VAR model 
%      p:               Number of lags in the VAR model                                                    (1 times 1)
%      horizons:        Number of horizons for the Impulse Response Functions (IRFs)                       (1 times 1)  
%                       (does not include the impact horizon 0)    
%      confidence:      Value for the standard and weak-IV robust confidence set                           (1 times 1)
%      figureorder:     Number representing the index of the next figure                                   (1 times 1)
%      Plugin:          Structure containing standard plug-in inference
%      InferenceMSW:    InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
%      time:            Time unit for the dataset (e.g. year, month, etc.)                                 (String)
%      columnnames:     Vector with the names for the endogenous variables, in the same order as ydata     (1 times n)
%      savdir:          Directory where the figures generated will be saved                                (String)
%      direct:          Directory where TestScript.m is located                                            (String) 
%      dataset_name:    The name of the dataset used for generating the figures (used in the output label) (String)
%      IRFselect:       Indices for the variables that the user wants separate IRF plots for
%      cumselect:       Indices for the variables that the user wants cumulative IRF plots for

% This version: Semptember 14th, 2018
% Last edited by José Luis Montiel-Olea

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

%check figureorder
if isempty(figureorder)
    
    error('figureorder must be assigned a value.');

end

if isnumeric(figureorder) == 0
    
    error('figureorder must be numeric.');

end

if numel(figureorder) ~= 1
    
    error('figureorder must have only one element');
    
end

if figureorder < 1
    
    error('figureorder must be >= 1');
    
end

%check Plugin.IRF

if isempty(Plugin.IRF)
    
    error('Plugin.IRF must not be empty.');

end

if isnumeric(Plugin.IRF) == 0
    
    error('Plugin.IRF must be numeric.');

end

if size(Plugin.IRF,1) ~= n
    
    error('Plugin.IRF must have n rows');
    
end

if size(Plugin.IRF,2) ~= horizons+1
    
    error('Plugin.IRF must have horizons + 1 columns');
    
end

%check Plugin.IRFstderror
if isempty(Plugin.IRFstderror)
    
    error('Plugin.IRFstderror must not be empty.');

end

if isnumeric(Plugin.IRFstderror) == 0
    
    error('Plugin.IRFstderror must be numeric.');

end

if size(Plugin.IRFstderror,1) ~= n
    
    error('Plugin.IRFstderror must have n rows');
    
end

if size(Plugin.IRFstderror,2) ~= horizons+1
    
    error('Plugin.IRFstderror must have horizons + 1 columns');
    
end

%check Plugin.IRFcum
if isempty(Plugin.IRFcum)
    
    error('Plugin.IRFcum must not be empty.');

end

if isnumeric(Plugin.IRFcum) == 0
    
    error('Plugin.IRFcum must be numeric.');

end

if size(Plugin.IRFcum,1) ~= n
    
    error('Plugin.IRFcum must have n rows');
    
end

if size(Plugin.IRFcum,2) ~= horizons+1
    
    error('Plugin.IRFcum must have horizons + 1 columns');
    
end

%check Plugin.IRFstderrorcum
if isempty(Plugin.IRFstderrorcum)
    
    error('Plugin.IRFstderrorcum must not be empty.');

end

if isnumeric(Plugin.IRFstderrorcum) == 0
    
    error('Plugin.IRFstderrorcum must be numeric.');

end

if size(Plugin.IRFstderrorcum,1) ~= n
    
    error('Plugin.IRFstderrorcum must have n rows');
    
end

if size(Plugin.IRFstderrorcum,2) ~= horizons+1
    
    error('Plugin.IRFstderrorcum must have horizons + 1 columns');
    
end

%check InferenceMSW.bootsIRFs

if isempty(InferenceMSW.bootsIRFs)
    
    error('InferenceMSW.bootsIRFs must not be empty.');

end

if isnumeric(InferenceMSW.bootsIRFs) == 0
    
    error('InferenceMSW.bootsIRFs must be numeric.');

end

if ndims(InferenceMSW.bootsIRFs) ~= 4
    
    error('InferenceMSW.bootsIRFs must be four dimensional');
    
end

if size(InferenceMSW.bootsIRFs,1) ~= n
    
    error('InferenceMSW.bootsIRFs''s first dimension must be n long');
    
end

if size(InferenceMSW.bootsIRFs,2) ~= horizons + 1
   
    error('InferenceMSW.bootsIRFs''s second dimension must be horizons + 1 long'); 
    
end

if size(InferenceMSW.bootsIRFs,3) ~= 2
   
    error('InferenceMSW.bootsIRFs''s third dimension must have size = 2'); 
    
end

if size(InferenceMSW.bootsIRFs,4) ~= 2
   
    error('InferenceMSW.bootsIRFs''s fourth dimension must have size = 2'); 
    
end

%check time
if ischar(time) == 0
    
    error('time must be a character array');
    
end

if isvector(time) ~= 1
    
    error('time must be a character vector.');
    
end

if size(time,1) ~= 1
    
    error('time must be a character row vector.');
    
end

%check columnnames
if isempty(columnnames)
    
    error('columnnames must be assigned a character array.');

end

if iscellstr(columnnames) == 0
    
    error('columnnames must be a cell array of character vectors.');

end

if size(columnnames,2) > n
    
     error('columnames size must be smaller or equal to the total amount of variables');
     
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

%check direct
if isempty(direct)
    
    direct('direct must be assigned a value.');

end

if ischar(direct) == 0
    
    error('direct must be a character array.');

end

if size(direct,1) > 1
    
    error('direct must only contain one string');
    
end

%check dataset_name
if ischar(dataset_name) == 0
    
    error('dataset_name must be a character array');
    
end

if isvector(dataset_name) ~= 1
    
    error('dataset_name must be a character vector.');
    
end

if size(time,1) ~= 1
    
    error('dataset_name must be a character row vector.');
    
end

%Check IRFselect
if isempty(IRFselect) ~= 0
   
    %IRFselect can be empty. If it isn't, then we need to check if it
    %follows the following criteria:
    if isnumeric(IRFselect) == 0
        
        error('IRFselect must be a numeric array');
        
    end
    
    if isvector(IRFselect) ~= 1
        
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
    
    if isvector(cumselect) ~= 1
        
        error('cumselect must be a row vector');
    
    end
    
    if size(cumselect,1) > 1
        
        error('cumselect must be a row vector');
        
    end
    
    if size(cumselect,2) > size(ydata, 2)
        
        error('cumselect size must be smaller or equal to the total amount of variables');
        
    end
    
end

end