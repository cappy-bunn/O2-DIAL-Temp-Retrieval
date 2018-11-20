filename = 'air.2018.nc';

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

air_name = 'air';
air_index = find(strcmp(var_names,air_name));
air_data = var_data{air_index};
air_dims = var_dims{air_index};

imagesc(squeeze(air_data(:,:,1)))

level_name = 'level';
level_index = find(strcmp(var_names,level_name));
level_data = var_data{level_index};
level_dims = var_dims{level_index};