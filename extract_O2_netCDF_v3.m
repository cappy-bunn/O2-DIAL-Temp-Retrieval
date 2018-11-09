

% Set maximum range [m] to use:
r_max = 6000;

% Set time to average over [mins]:
movavg_t = 20;

% Set range to average over [m]:
movavg_r = 150;


% ======== Import NetCDF data ========
ncid = netcdf.open('wv_dial05.181020.Python.nc');

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


% ----- Time and Range profiles -----

% Temperature Time and Range profiles
t_name = dim_names{30};
t_ind = find(strcmp(var_names, t_name));
t_T = double(var_data{t_ind});                  % Time [sec]
t_hr_T = t_T./3600;                             % Time [hr]

r_name = dim_names{31};
r_index = find(strcmp(var_names, r_name));
r_T = double(var_data{r_index});
ind_r_T = r_T <= r_max;
r_T = r_T(ind_r_T);                             % Range up to r_max [m]

% Pressure Time and Range profiles
t_name = dim_names{32};
t_ind = find(strcmp(var_names, t_name));
t_P = double(var_data{t_ind});                  % Time [sec]
t_hr_P = t_P./3600;                             % Time [hr]

r_name = dim_names{33};
r_index = find(strcmp(var_names, r_name));
r_P = double(var_data{r_index});
ind_r_P = r_P <= r_max;
r_P = r_P(ind_r_P);                             % Range up to r_max [m]

% Online O2 Time and Range profiles
t_name = dim_names{19};
t_ind = find(strcmp(var_names, t_name));
t_o2on = double(var_data{t_ind});%#ok<*FNDSB>   % Time [sec] 
t_hr_o2on = t_o2on./3600;                       % Time [hr]

r_name = dim_names{20};
r_index = find(strcmp(var_names, r_name));
r_o2on = double(var_data{r_index});
ind_r_o2on = r_o2on >= r_T(2) & r_o2on<= r_max;
r_o2on = r_o2on(ind_r_o2on);                    % Range from 2nd element of r_T up to r_max 

% Offline O2 Time and Range profiles (same as online)
t_name = dim_names{22};
t_ind = find(strcmp(var_names, t_name));
t_o2off = double(var_data{t_ind});%#ok<*FNDSB>  % Time [sec] 
t_hr_o2off = t_o2off./3600;                     % Time [hr]

r_name = dim_names{23};
r_index = find(strcmp(var_names, r_name));
r_o2off = double(var_data{r_index});
ind_r_o2off = r_o2off >= r_T(2) & r_o2off<= r_max;
r_o2off = r_o2off(ind_r_o2off);                 % Range from 2nd element of r_T up to r_max

% Absolute Humidity Time and Range profiles
t_name = dim_names{24};
t_ind = find(strcmp(var_names, t_name));
t_wv = double(var_data{t_ind});                 % Time [sec]
t_hr_wv = t_wv./3600;                           % Time [hr]

r_name = dim_names{25};
r_index = find(strcmp(var_names, r_name));
r_wv = double(var_data{r_index});
ind_r_wv = r_wv >= r_T(2) & r_wv<= r_max;
r_wv = r_wv(ind_r_wv);                          % Range from 2nd element of r_T up to r_max

% Denoised Aerosol Backscatter Coefficient Time and Range profiles
t_name = dim_names{38};
t_ind = find(strcmp(var_names, t_name));
t_ABC = double(var_data{t_ind});                % Time [sec]
t_hr_ABC = t_ABC./3600;                         % Time [hr]

r_name = dim_names{39};
r_index = find(strcmp(var_names, r_name));
r_ABC = double(var_data{r_index});
ind_r_ABC = r_ABC <= r_max;
r_ABC = r_ABC(ind_r_ABC);                       % Range up to r_max [m]

% Denoised Backscatter Ratio Time and Range profiles
t_name = dim_names{42};
t_ind = find(strcmp(var_names, t_name));
t_BR = double(var_data{t_ind});                 % Time [sec]
t_hr_BR = t_BR./3600;                           % Time [hr]

r_name = dim_names{43};
r_index = find(strcmp(var_names, r_name));
r_BR = double(var_data{r_index});
ind_r_BR = r_BR <= r_max;
r_BR = r_BR(ind_r_BR);                          % Range up to r_max [m]


% ----- Read in Temperature and Pressure profiles -----

T_name = 'Temperature';
T_index = find(strcmp(var_names, T_name));
T_data = var_data{T_index};
T_dims = var_dims{T_index};
T = double(T_data);
T = T(ind_r_T,:);                               % Temperature [K]

P_name = 'Pressure';
P_index = find(strcmp(var_names, P_name));
P_data = var_data{P_index};
P_dims = var_dims{P_index};
P = double(P_data);                             
P = P(ind_r_P,:);                               % Pressure [atm]

% ----- O2 Backscatter Channel Raw Data -----------------

% ONLINE
o2on_name = 'O2_Online_Backscatter_Channel_Raw_Data';
o2on_index = find(strcmp(var_names, o2on_name));
o2on_data = var_data{o2on_index};
o2on_dims = var_dims{o2on_index};
o2on = double(o2on_data);
o2on = o2on(ind_r_o2on,:);                      % Online return signal

% Moving average in range. Range bin is 37.5 meters.
window_r = round(movavg_r/37.5);                % Number of elements to average over
o2on = movmean(o2on,window_r,1,'Endpoints','discard');

% Moving average in time. Data is already averaged over 5 minutes.
window_t = round(movavg_t/5);                   % Number of elements to average over
o2on = movmean(o2on,window_t,2,'Endpoints','discard');


% OFFLINE
o2off_name = 'O2_Offline_Backscatter_Channel_Raw_Data';
o2off_index = find(strcmp(var_names, o2off_name));
o2off_data = var_data{o2off_index};
o2off_dims = var_dims{o2off_index};
o2off = double(o2off_data);
o2off = o2off(ind_r_o2off,:);                   % Offline return signal

% ----- Water Vapor Number Density --------------------------

wv_name = 'Absolute_Humidity';
wv_index = find(strcmp(var_names, wv_name));
wv_data = var_data{wv_index};
wv_dims = var_dims{wv_index};

wvm_name = 'Absolute_Humidity_mask';
wvm_index = find(strcmp(var_names, wvm_name));
wvm_data = var_data{wvm_index};
wvm_dims = var_dims{wvm_index};
wvm_data = double(wvm_data);

wv_profile = wv_data.*abs((1-wvm_data));    % Apply mask
wv_profile(wv_profile<0) = 0;
wv_profile = double(wv_profile);
wv_profile = wv_profile(ind_r_wv,:);        % WV profile [g/m^3]
m_h2o = 18.015/6.022E23;                    % Mass of H2O molecule [g]

wv_num = wv_profile./m_h2o;                 % WV number density [#/m^3]


% ----- Denoised Aerosol Backscatter Coefficient ---------------------

dabc_name = 'Denoised_Aerosol_Backscatter_Coefficient';
dabc_index = find(strcmp(var_names, dabc_name));
dabc_data = var_data{dabc_index};
dabc_dims = var_dims{dabc_index};

dabcm_name = 'Denoised_Aerosol_Backscatter_Coefficient_mask';
dabcm_index = find(strcmp(var_names, dabcm_name));
dabcm_data = var_data{dabcm_index};
dabcm_dims = var_dims{dabcm_index};
dabcm_data = double(dabcm_data);

betaa = double(dabc_data.*(1-dabcm_data));     % Apply mask
betaa = betaa(ind_r_ABC,:);                    % Aerosol Backscatter Coefficient


% ----- Denoised Backscatter Ratio ----------------------------

dbr_name = 'Denoised_Backscatter_Ratio';
dbr_index = find(strcmp(var_names, dbr_name));
dbr_data = var_data{dbr_index};
dbr_dims = var_dims{dbr_index};

dbrm_name = 'Denoised_Aerosol_Backscatter_Coefficient_mask';
dbrm_index = find(strcmp(var_names, dbrm_name));
dbrm_data = var_data{dbrm_index};
dbrm_dims = var_dims{dbrm_index};
dbrm_data = double(dbrm_data);

br = double(dbr_data.*(1-dbrm_data));        % Apply mask
br = br(ind_r_BR,:);                         % Backscatter Ratio

