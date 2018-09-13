function SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name)

if ~strncmp(version, '9', 1)
    
    warning('You are using a version other than v9, in which this program was tested.')

end 

%check p value
if isempty(p)
    
    error('p should be assigned a value.');

end

if isnumeric(p) == 0
    
    error('p should be numeric.');

end

if numel(p) ~= 1
    
    error('p should have only one element');
    
end

if floor(p) ~= p
    
    error('p should be an integer.');
    
end

%check confidence value

if isempty(confidence)
    
    error('confidence should be assigned a value.');

end

if isnumeric(confidence) == 0
    
    error('confidence should be numeric.');

end

if numel(confidence) ~= 1
    
    error('confidence should have only one element');
    
end

if(confidence <= 0 || confidence >= 1)
    
    error('confidence should be > 0 and < 1');
    
end

%check ydata
if isempty(ydata) == 1
    
    error('ydata should be a (T times n) matrix. It is now empty');

end

if isnumeric(ydata) == 0
    
    error('ydata should be numeric.');
    
end

if ismatrix(ydata) == 0
    
    error('ydata should be two dimensional (T times n).');

end 

if size(ydata,1) < size(ydata,2)
    
    warning('ydata: number of rows < number of columns. ydata should be a T times n.');
    
end

T = size(ydata,1);

%check z
if isempty(z) == 1
    
    error('z should be (T times 1) vector. It is now empty');

end

if isnumeric(z) == 0
    
    error('z should be numeric');
    
end

if ismatrix(z) == 0
    
    error('z should be two dimensional (T times 1)');
    
end

if size(z,2) ~= 1
    
    error('There should be only one instrument (z should have only one column) should be');
    
end

if size(z,1) ~= T
    
    error('z should have T rows');
    
end

%check NWlags
if isempty(NWlags)
    
    error('NWlags should be assigned a value.');

end

if isnumeric(NWlags) == 0
    
    error('NWlags should be numeric.');

end

if numel(NWlags) ~= 1
    
    error('NWlags should have only one element');
    
end

if floor(NWlags) ~= NWlags
    
    error('NWlags should be an integer.');
    
end

%Check norm
if isempty(norm)
    
    error('norm should be assigned a value.');

end

if isnumeric(norm) == 0
    
    error('norm should be numeric.');

end

if numel(norm) ~= 1
    
    error('norm should have only one element');
    
end

if floor(norm) ~= norm
    
    error('norm should be an integer.');
    
end

if norm > size(ydata, 2)
    
    error('norm should be smaller or equal to the total amount of variables');
    
end

if norm < 1
    
    error('norm should be >= 1');
    
end

%check scale
if isempty(scale)
    
    error('scale should be assigned a value.');

end

if isnumeric(scale) == 0
    
    error('scale should be numeric.');

end

if numel(scale) ~= 1
    
    error('scale should have only one element');
    
end

if scale <= 0
    
    error('scale should be > 0.');

end

if floor(scale) ~= scale
    
    error('scale should be an integer.');
    
end
    
%check horizons

if isempty(horizons)
    
    error('horizons should be assigned a value.');

end

if isnumeric(horizons) == 0
    
    error('horizons should be numeric.');

end

if numel(horizons) ~= 1
    
    error('horizons should have only one element');
    
end

if horizons < 0
    
    error('horizons should be >= 0.');
    
end

if floor(horizons) ~= horizons
    
    error('horizons should be an integer.');
    
end

%check savdir
if isempty(savdir)
    
    savdir('savdir should be assigned a value.');

end

if ischar(savdir) == 0
    
    error('savdir should be a character array.');

end

if size(savdir,1) > 1
    
    error('savdir should only contain one string');
    
end

%columnnames
if isempty(columnnames)
    
    error('columnnames should be assigned a character array.');

end

if iscellstr(columnnames) == 0
    
    error('columnnames should be a cell array of character vectors.');

end

if size(columnnames,2) > size(ydata, 2)
    
     error('columnames size should be smaller or equal to the total amount of variables');
     
end

%Check IRFselect
if isempty(IRFselect) ~= 0
   
    %IRFselect can be empty. If it isn't, then we need to check if it
    %follows the following criteria:
    if isnumeric(IRFselect) == 0
        
        error('IRFselect should be a numeric array');
        
    end
    
    if ismatrix(IRFselect) ~= 2
        
        error('IRFselect should be a row vector');
    
    end
    
    if size(IRFselect,1) > 1
        
        error('IRFselect should be a row vector');
        
    end
    
    if size(IRFselect,2) > size(ydata, 2)
        
        error('IRFselect size should be smaller or equal to the total amount of variables');
        
    end
    
    
    
end

%Check cumselect
if isempty(cumselect) ~= 0
   
    %cumselect can be empty. If it isn't, then we need to check if it
    %follows the following criteria:
    if isnumeric(cumselect) == 0
        
        error('cumselect should be a numeric array');
        
    end
    
    if ismatrix(cumselect) ~= 2
        
        error('cumselect should be a row vector');
    
    end
    
    if size(cumselect,1) > 1
        
        error('cumselect should be a row vector');
        
    end
    
    if size(cumselect,2) > size(ydata, 2)
        
        error('cumselect size should be smaller or equal to the total amount of variables');
        
    end
    
end
    
%check time

if ischar(time) == 0
    
    error('time should be a character array');
    
end

if ismatrix(time) == 0
    
    error('time should be a character array. Dimension is > 2');
    
end

if size(time, 1) > 1
    
    error('time should be a character array. It is a matrix');
    
end

%check dataset_name

    
if ischar(dataset_name) == 0
    
    error('dataset_name should be a character array');
    
end

if ismatrix(dataset_name) == 0
    
    error('dataset_name should be a character array. Dimension is > 2');
    
end

if size(dataset_name, 1) > 1
    
    error('dataset_name should be a character array. It is a matrix');
    
end

end

    



    


 
    
    





