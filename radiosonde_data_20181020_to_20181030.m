filename = 'raob_soundings24581.cdf';

% Open netcdf file
ncid = netcdf.open(filename);

% Read the number of dimensions and number of variables in the file
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

% Read in information about each dimension
dim_names = {};
dim_sz = [];
data = [];
for i = 0:numdims-1

    [dimname, dimlen] = netcdf.inqDim(ncid,i);
    dim_names{i+1} = dimname;
    dim_sz = [dim_sz, dimlen];
    
end

% Read variables
var_names = {};
var_types = [];
var_dims = {};
var_natts = [];
var_data = {};
for i =0:numvars-1
   
    % Inquire about variables
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i);  
    var_names{i+1} = varname;
    var_types = [var_types, xtype];
    var_dims{i+1} = dimids;
    var_natts = [var_natts, natts];
    
    % Read in data
    count = dim_sz(dimids+1);
    start = 0 * count;
    data = netcdf.getVar(ncid,i,start,count);
    var_data{i+1} = data;
        
end

% Time/Date
time_name = 'synTime';
time_index = find(strcmp(var_names,time_name));
time_data = var_data{time_index};
time_dims = var_dims{time_index};

time_data = datetime(time_data,'convertfrom','posixtime');  % Convert from unix time to normal time/date

% Height [m]
height_name = 'htMan';
height_index = find(strcmp(var_names,height_name));
height_data = var_data{height_index};
height_dims = var_dims{height_index};

height_data(height_data == max(max(height_data))) = NaN;
height_data(1,:) = 0;   % Change first row from station elevation to 0 (surface)

% Mandatory levels
tempMan_name = 'tpMan';
tempMan_index = find(strcmp(var_names,tempMan_name));
tempMan_data = var_data{tempMan_index};
tempMan_dims = var_dims{tempMan_index};

tempMan_data(tempMan_data == max(max(tempMan_data))) = NaN;
tempMan_data(tempMan_data == max(max(tempMan_data))) = NaN;

% Significant levels wrt Temperature
tempST_name = 'tpSigT';
tempST_index = find(strcmp(var_names,tempST_name));
tempST_data = var_data{tempST_index};
tempST_dims = var_dims{tempST_index};

tempST_data(tempST_data == max(max(tempST_data))) = NaN;

plot(tempMan_data(1:5,12),height_data(1:5,12))

tempST_name = 'sigTLevel';
tempST_index = find(strcmp(dim_names,tempST_name));
tempST_data = var_data{tempST_index};
tempST_dims = var_dims{tempST_index};