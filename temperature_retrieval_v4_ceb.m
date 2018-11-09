% ========================================================================
% Program for retrieving temperature profiles from initial O2 data
% ========================================================================
close all;
clear;
%clc;

%% ========================================================================
% Input parameters for the netcdf file retrieval
% =========================================================================
filename = 'wv_dial05.181030.Python.nc';
r_max = 6000;           % Maximum range [m]
movavg_r = 300;         % Set range to average over [m] should be a multiple of 37.5 m
movavg_t = 20;          % Set time to average over [mins] should be a multiple of 5 minutes
nighttime = 1;          % 1 = YES, only look at nighttime; 0 = NO, look at full day
dim_offset = 1;         % Set offset in 'dim_names'

% ========================================================================
% Reading in data from netcdf files
% ========================================================================
% ======== Import NetCDF data ========
% ========================================================================

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


% ========== Data strings differ starting on Oct 25th, 2018 ===============

% Search filename for the date
date = str2double(regexp(filename,'\d{6}','match','once'));

if date > 181024
    % Data 'Temperature' and 'Pressure' not available. Use these instead:
    T_str = 'Temperature_HSRL';
    P_str = 'Pressure_HSRL';
    
    % --- Surface Temperature and Pressure ---
    % Surface Temperature is split off from the rest of the temperature profile
    T_surf_name = 'Surface_Temperature_HSRL';
    T_surf_index = find(strcmp(var_names,T_surf_name));
    T_surf_data = var_data{T_surf_index};
    T_surf_dims = var_dims{T_surf_index};
    T_surf = double(T_surf_data);                   % Surface Temperature [K]
    
    % Surface temperature time and range profiles
    t_name = dim_names{T_surf_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names, t_name));
    t_T_surf = double(var_data{t_ind});             % Time [sec]
    
    r_name = dim_names{T_surf_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names, r_name));
    r_T_surf = double(var_data{r_ind});             % Range = 0 m

    % Surface Pressure is split off from the rest of the pressure profile
    P_surf_name = 'Surface_Pressure_HSRL';
    P_surf_index = find(strcmp(var_names,P_surf_name));
    P_surf_data = var_data{P_surf_index};
    P_surf_dims = var_dims{P_surf_index};
    P_surf = double(P_surf_data);                   % Surface Pressure [atm]
    
    % Surface Pressure time and range profiles
    t_name = dim_names{P_surf_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names, t_name));
    t_P_surf = double(var_data{t_ind});             % Time [sec]
    
    r_name = dim_names{P_surf_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names, r_name));
    r_P_surf = double(var_data{r_ind});             % Range = 0 m
    
    
    % --- Temperature and Pressure from DIAL 4 ---
    % Extract 'Temperature' and 'Pressure' from another DIAL to compare
    
    filename2 = strrep(filename,'wv_dial05','wv_dial04');
    
    % Open netcdf file
    ncid2 = netcdf.open(filename2);

    % Read the number of dimensions and number of variables in the file
    [numdims2, numvars2, numglobalatts2, unlimdimID2] = netcdf.inq(ncid2);

    % Read in information about each dimension
    dim_names2 = {};
    dim_sz2 = [];
    data2 = [];
    for i = 0:numdims2-1

        [dimname2, dimlen2] = netcdf.inqDim(ncid2,i);
        dim_names2{i+1} = dimname2;
        dim_sz2 = [dim_sz2, dimlen2];
    
    end

    % Read variables
    var_names2 = {};
    var_types2 = [];
    var_dims2 = {};
    var_natts2 = [];
    var_data2 = {};
    
    for i =0:numvars2-1
   
        % Inquire about variables
        [varname2,xtype2,dimids2,natts2] = netcdf.inqVar(ncid2,i);  
        var_names2{i+1} = varname2;
        var_types2 = [var_types2, xtype2];
        var_dims2{i+1} = dimids2;
        var_natts2 = [var_natts2, natts2];
    
        % Read in data
        count2 = dim_sz2(dimids2+1);
        start2 = 0 * count2;
        data2 = netcdf.getVar(ncid2,i,start2,count2);
        var_data2{i+1} = data2;
        
    end
    
    T2_name = 'Temperature';
    T2_index = find(strcmp(var_names2, T2_name));
    T2_data = var_data2{T2_index};
    T2_dims = var_dims2{T2_index};
    T2 = double(T2_data);

    % Temperature time and range profiles
    t_name = dim_names2{T2_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names2, t_name));
    t_T2 = double(var_data2{t_ind});                   % Time [sec]

    r_name = dim_names2{T2_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names2, r_name));
    r_T2 = double(var_data2{r_ind});
    ind_r_T2 = r_T2 <= r_max;                         
    r_T2 = r_T2(ind_r_T2);                             % Range up to r_max [m]

    T2 = T2(ind_r_T2,:);                               % Temperature [K]

    P2_name = 'Pressure';
    P2_index = find(strcmp(var_names2, P2_name));
    P2_data = var_data2{P2_index};
    P2_dims = var_dims2{P2_index};
    P2 = double(P2_data);                             

    % Pressure time and range profiles
    t_name = dim_names2{P2_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names2, t_name));
    t_P2 = double(var_data2{t_ind});                  % Time [sec]

    r_name = dim_names2{P2_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names2, r_name));
    r_P2 = double(var_data2{r_ind});
    ind_r_P2 = r_P2 <= r_max;
    r_P2 = r_P2(ind_r_P2);                             % Range up to r_max [m]

    P2 = P2(ind_r_P2,:);                               % Pressure [atm]
else
    T_str = 'Temperature';
    P_str = 'Pressure';
end


% ================= Temperature and Pressure profiles =====================

T_name = T_str;
T_index = find(strcmp(var_names, T_name));
T_data = var_data{T_index};
T_dims = var_dims{T_index};
T = double(T_data);

% Temperature time and range profiles
t_name = dim_names{T_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_T = double(var_data{t_ind});                  % Time [sec]

r_name = dim_names{T_dims(1)+dim_offset};
r_ind = find(strcmp(var_names, r_name));
r_T = double(var_data{r_ind});
ind_r_T = r_T > 0 & r_T <= r_max;                         
r_T = r_T(ind_r_T);                             % Range up to r_max [m] (not including surface)

if date <= 181024   % Split off surface if not already done.
    T_surf = T(1,:);                            % Surface Temperature [K]
end

T = T(ind_r_T,:);                               % Temperature [K] (not including surface T)

P_name = P_str;
P_index = find(strcmp(var_names, P_name));
P_data = var_data{P_index};
P_dims = var_dims{P_index};
P = double(P_data);                             

% Pressure time and range profiles
t_name = dim_names{P_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_P = double(var_data{t_ind});                  % Time [sec]

r_name = dim_names{P_dims(1)+dim_offset};
r_ind = find(strcmp(var_names, r_name));
r_P = double(var_data{r_ind});
ind_r_P = r_P > 0 & r_P <= r_max;
r_P = r_P(ind_r_P);                             % Range up to r_max [m] (not including surface)

if date <= 181024   % Split off surface if not already done.
    P_surf = P(1,:);                            % Surface Pressure [atm]
end

P = P(ind_r_P,:);                               % Pressure [atm] (not including surface P)


% ================= O2 Backscatter Channel ================================

% ----- O2 ONLINE ---------------------------------------------------------

o2on_name = 'O2_Online_Backscatter_Channel_Raw_Data';
o2on_index = find(strcmp(var_names, o2on_name));
o2on_data = var_data{o2on_index};
o2on_dims = var_dims{o2on_index};
o2on = double(o2on_data);

% --- Online O2 time and range profiles ---
t_name = dim_names{o2on_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_o2on = double(var_data{t_ind});%#ok<*FNDSB>   % Time [sec] 

r_name = dim_names{o2on_dims(1)+dim_offset};
r_index = find(strcmp(var_names, r_name));
r_o2on_prep = double(var_data{r_index});
ind_r_o2on_var = r_o2on_prep >= r_T(1) & r_o2on_prep <= r_max;
r_o2on = r_o2on_prep(ind_r_o2on_var);           % Range [m] up to r_max 

% Set max r to 15 km so background counts may be used
ind_r_o2on_15 = r_o2on_prep > 0 & r_o2on_prep <= 15000;
o2on = o2on(ind_r_o2on_15,:);                   % Online return signal to 15 km


% --- Moving average in RANGE ---
window_r = round(movavg_r/37.5);                % Number of elements to average over
row1_o2on = o2on(1,:);                          % First row of O2 online backscatter
endrow_o2on = o2on(end,:);                      % Last row of O2 online backscatter
o2on_prep = o2on;

% Appending the correct number of first/last rows depending on if the RANGE window is even or odd.
if mod(window_r,2) == 0
    % Concatenate repeated beginning row(s) if even
    for i = 1:((window_r/2)-1)
        o2on_prep = [row1_o2on; o2on_prep];
    end
    % Concatenate repeated end row(s) if even
    for i = 1:(window_r/2)
        o2on_prep = [o2on_prep; endrow_o2on]; 
    end
else
    % Concatenate repeated beginning & end row(s) if odd
    for i = 1:((window_r-1)/2)
        o2on_prep = [row1_o2on; o2on_prep; endrow_o2on];
    end
end

% Online O2 averaged over RANGE
o2on_ravg = movmean(o2on_prep,window_r,1,'Endpoints','discard');    


% --- Moving average in TIME ---
window_t = round(movavg_t/5);                   % Number of elements to average over
col1_o2on = o2on(:,1);                          % First column of O2 online backscatter
endcol_o2on = o2on(:,end);                      % Last column of O2 online backscatter
o2on_prep = o2on_ravg;

% Appending the correct number of first/last columns depending on if the TIME window is even or odd.
if mod(window_t,2) == 0
    % Concatenate repeated beginning columns(s) if even
    for i = 1:((window_t/2)-1)
        o2on_prep = [col1_o2on o2on_prep];
    end
    % Concatenate repeated end columns(s) if even
    for i = 1:(window_t/2)
        o2on_prep = [o2on_prep endcol_o2on]; 
    end
else
    % Concatenate repeated beginning & end columns(s) if odd
    for i = 1:((window_t-1)/2)
        o2on_prep = [col1_o2on o2on_prep endcol_o2on];
    end
end

% Online O2 averaged over time and range
o2on_rtavg = movmean(o2on_prep,window_t,2,'Endpoints','discard');
                      

% --- Subtract background (last 20 bins) ---
background_o2on = repmat(mean(o2on_rtavg(end-20:end,:),1),size(o2on_rtavg,1),1);
o2on_avg_15 = o2on_rtavg - background_o2on; % Background subtracted, averaged O2 Online counts (range 15 km)

% --- Adjusts range from 15 km back to r_max ---
o2on_avg = o2on_avg_15(ind_r_o2on_var,:);      

% --- Interpolate if necessary ---
% If averaging windows are even, the O2 counts were evaluated half a bin up in range and time
% Interpolate counts to place them back at original locations
if mod(window_r,2) == 0
    r_o2on_shift = r_o2on + (mean(diff(r_o2on))/2);
    o2on_avg = interp1(r_o2on_shift,o2on_avg,r_o2on,'linear','extrap');
end

if mod(window_t,2) == 0
    t_o2on_shift = t_o2on + (mean(diff(t_o2on))/2);
    o2on_avg_transpose = interp1(t_o2on_shift,o2on_avg',t_o2on,'linear','extrap');
    o2on_avg = o2on_avg_transpose';
end

o2on_avg(o2on_avg <= 0) = NaN;

% ----- O2 OFFLINE --------------------------------------------------------

o2off_name = 'O2_Offline_Backscatter_Channel_Raw_Data';
o2off_index = find(strcmp(var_names, o2off_name));
o2off_data = var_data{o2off_index};
o2off_dims = var_dims{o2off_index};
o2off = double(o2off_data);

% --- Offline O2 time and range profiles ---
t_name = dim_names{o2off_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_o2off = double(var_data{t_ind});%#ok<*FNDSB>  % Time [sec] 

r_name = dim_names{o2off_dims(1)+dim_offset};
r_ind = find(strcmp(var_names, r_name));
r_o2off_prep = double(var_data{r_ind});
ind_r_o2off_var = r_o2off_prep >= r_T(1) & r_o2off_prep <= r_max;
r_o2off = r_o2off_prep(ind_r_o2off_var);        % Range [m] up to r_max

% Set max r to 15 km so background counts may be used
ind_r_o2off_15 = r_o2off_prep > 0 & r_o2off_prep <= 15000; 
o2off = o2off(ind_r_o2off_15,:);                % Offline return signal to 15 km


% --- Moving average in RANGE ---
row1_o2off = o2off(1,:);                        % First row of O2 offline backscatter
endrow_o2off = o2off(end,:);                    % Last row of O2 offline backscatter
o2off_prep = o2off;

% Appending the correct number of first/last rows depending on if the RANGE window is even or odd.
if mod(window_r,2) == 0
    % Concatenate repeated beginning row(s) if even
    for i = 1:((window_r/2)-1)
        o2off_prep = [row1_o2off; o2off_prep];
    end
    % Concatenate repeated end row(s) if even
    for i = 1:(window_r/2)
        o2off_prep = [o2off_prep; endrow_o2off]; 
    end
else
    % Concatenate repeated beginning & end row(s) if odd
    for i = 1:((window_r-1)/2)
        o2off_prep = [row1_o2off; o2off_prep; endrow_o2off];
    end
end

% Offline O2 averaged over RANGE
o2off_ravg = movmean(o2off_prep,window_r,1,'Endpoints','discard');    


% --- Moving average in TIME ---
col1_o2off = o2off(:,1);           % First column of O2 offline backscatter
endcol_o2off = o2off(:,end);       % Last column of O2 offline backscatter
o2off_prep = o2off_ravg;

% Appending the correct number of first/last columns depending on if the TIME window is even or odd.
if mod(window_t,2) == 0
    % Concatenate repeated beginning columns(s) if even
    for i = 1:((window_t/2)-1)
        o2off_prep = [col1_o2off o2off_prep];
    end
    % Concatenate repeated end columns(s) if even
    for i = 1:(window_t/2)
        o2off_prep = [o2off_prep endcol_o2off]; 
    end
else
    % Concatenate repeated beginning & end columns(s) if odd
    for i = 1:((window_t-1)/2)
        o2off_prep = [col1_o2off o2off_prep endcol_o2off];
    end
end

% Offline O2 averaged over time and range
o2off_rtavg = movmean(o2off_prep,window_t,2,'Endpoints','discard');


% --- Subtract background (average of last 20 bins) ---
background_o2off = repmat(mean(o2off_rtavg(end-20:end,:),1),size(o2off_rtavg,1),1);

o2off_avg_15 = o2off_rtavg - background_o2off; % Background subtracted, averaged O2 Offline counts (15 km)

% --- Adjust range from 15 km back to r_max ---
o2off_avg = o2off_avg_15(ind_r_o2off_var,:);

% --- Interpolate if necessary ---
% If averaging windows are even, the O2 counts were evaluated half a bin up in range and time
% Interpolate counts to place them back at original locations
if mod(window_r,2) == 0
    r_o2off_shift = r_o2off + (mean(diff(r_o2off))/2);
    o2off_avg = interp1(r_o2off_shift,o2off_avg,r_o2off,'linear','extrap');
end

if mod(window_t,2) == 0
    t_o2off_shift = t_o2off + (mean(diff(t_o2off))/2);
    o2off_avg_transpose = interp1(t_o2off_shift,o2off_avg',t_o2off,'linear','extrap');
    o2off_avg = o2off_avg_transpose';
end

%       ========================================================
%       ==== Final O2 counts are 'o2on_avg' and 'o2off_avg' ====
%       ========================================================


% ---- O2 PLOTS ----
% Plot original O2 online and background subtracted, averaged O2 online
figure
imagesc(o2on(ind_r_o2on_var,:))
set(gca,'YDir','normal')
figure
imagesc(o2on_avg)
set(gca,'YDir','normal')
                                    

% ============= Water Vapor Number Density ================================

% Must get 'Absolute Humidity' from another DIAL if date is after Oct 24th.
if date > 181024
    wv_name = 'Absolute_Humidity';
    wv_index = find(strcmp(var_names2, wv_name));
    wv_data = var_data2{wv_index};
    wv_dims = var_dims2{wv_index};

    wvm_name = 'Absolute_Humidity_mask';
    wvm_index = find(strcmp(var_names2, wvm_name));
    wvm_data = double(var_data2{wvm_index});

    % Absolute Humidity time and range profiles
    t_name = dim_names2{wv_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names2, t_name));
    t_wv = var_data2{t_ind};             

    r_name = dim_names2{wv_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names2, r_name));
    r_wv = double(var_data2{r_ind});
else
    wv_name = 'Absolute_Humidity';
    wv_index = find(strcmp(var_names, wv_name));
    wv_data = var_data{wv_index};
    wv_dims = var_dims{wv_index};

    wvm_name = 'Absolute_Humidity_mask';
    wvm_index = find(strcmp(var_names, wvm_name));
    wvm_data = double(var_data{wvm_index});

    % Absolute Humidity time and range profiles
    t_name = dim_names{wv_dims(2)+dim_offset};
    t_ind = find(strcmp(var_names, t_name));
    t_wv = var_data{t_ind};

    r_name = dim_names{wv_dims(1)+dim_offset};
    r_ind = find(strcmp(var_names, r_name));
    r_wv = double(var_data{r_ind});
end

wv_profile = wv_data.*((1-wvm_data));       % Apply mask
wv_profile(wv_profile<0) = 0;
wv_profile = double(wv_profile);

t_wv = double(t_wv);                        % Time [sec]
ind_r_wv = r_wv >= r_T(1) & r_wv<= r_max;
r_wv = r_wv(ind_r_wv);                      % Range up to r_max

wv_profile = wv_profile(ind_r_wv,:);        % WV profile [g/m^3]
m_h2o = 18.015/6.022E23;                    % Mass of H2O molecule [g]

wv_num = wv_profile./m_h2o;                 % WV number density [#/m^3]


% =========== Denoised Aerosol Backscatter Coefficient ====================

dabc_name = 'Denoised_Aerosol_Backscatter_Coefficient';
dabc_index = find(strcmp(var_names, dabc_name));
dabc_data = var_data{dabc_index};
dabc_dims = var_dims{dabc_index};

dabcm_name = 'Denoised_Aerosol_Backscatter_Coefficient_mask';
dabcm_index = find(strcmp(var_names, dabcm_name));
dabcm_data = double(var_data{dabcm_index});

% Denoised Aerosol Backscatter Coefficient time and range profiles
t_name = dim_names{dabc_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_ABC = double(var_data{t_ind});                % Time [sec]

r_name = dim_names{dabc_dims(1)+dim_offset};
r_ind = find(strcmp(var_names, r_name));
r_ABC = double(var_data{r_ind});
ind_r_ABC = r_ABC <= r_max;
r_ABC = r_ABC(ind_r_ABC);                       % Range up to r_max [m]

betaa = double(dabc_data.*(1-dabcm_data));      % Apply mask

betaa = betaa(ind_r_ABC,:);                     % Aerosol Backscatter Coefficient


% ============== Denoised Backscatter Ratio ===============================

dbr_name = 'Denoised_Backscatter_Ratio';
dbr_index = find(strcmp(var_names, dbr_name));
dbr_data = var_data{dbr_index};
dbr_dims = var_dims{dbr_index};

dbrm_name = 'Denoised_Backscatter_Ratio_mask';
dbrm_index = find(strcmp(var_names, dbrm_name));
dbrm_data = double(var_data{dbrm_index});

% Denoised Backscatter Ratio time and range profiles
t_name = dim_names{dbr_dims(2)+dim_offset};
t_ind = find(strcmp(var_names, t_name));
t_BR = double(var_data{t_ind});                 % Time [sec]
t_hr_BR = t_BR./3600;                           % Time [hr]

r_name = dim_names{dbr_dims(1)+dim_offset};
r_ind = find(strcmp(var_names, r_name));
r_BR = double(var_data{r_ind});
ind_r_BR = r_BR <= r_max;
r_BR = r_BR(ind_r_BR);                          % Range up to r_max [m]

br = double(dbr_data.*(1-dbrm_data));           % Apply mask

br = br(ind_r_BR,:);                            % Backscatter Ratio

% =========================================================================
% The variables from the above are as follows:
% T = temperature in K
% T_surf = surface temperature in K
% P = pressure in atm
% P_surf = surface pressure in atm
% T2 = Temperature from DIAL 4 if interested (surf T included)
% P2 = Pressure from DIAL 4 if interested (surf P included)
% o2on_avg = O2 online return counts averaged over range and time
% o2off_avg = O2 offline return counts averaged over range and time
% wv_num = water vapor number density
% betaa = denoised aerosol backscatter coefficient
% br = backscatter ratio
% =========================================================================
% Renaming variables for second part of program
rm = r_o2on;
rkm = rm./1000;
rangebin = movavg_r;
tsec = t_o2on;
thr = tsec./3600; 
%% ========================================================================
%  The following part of the program completes the temperature inversion
% =========================================================================
% =========================================================================
%  input parameters for the second portion of the program
if nighttime == 1
    maxtpro = 150;              % Use t profiles up to maxtpro (afterwards is daytime)
    tsec = tsec(1:maxtpro);
    thr = thr(1:maxtpro);
    T = T(:,1:maxtpro);
    T_surf = T_surf(1:maxtpro);
    P = P(:,1:maxtpro);
    P_surf = P_surf(1:maxtpro);
    o2on_avg = o2on_avg(:,1:maxtpro);
    o2off_avg = o2off_avg(:,1:maxtpro);
    wv_num = wv_num(:,1:maxtpro);
    betaa = betaa(:,1:maxtpro);
    br = br(:,1:maxtpro);
    
end


c = 3.00e8;                 % speed of light in m/s
kb = 1.38065e-23;           % Boltzman's constant in J/K
h = 6.626e-34;              % Planck's constant in Js
No = 2.687e25;              % Loschmidt's number in 1/m^3 referenced to 273.15 k and 1 atm
mo2 = 5.314e-26;            % mass o2 molecule in kg

l_o2_on = 769.795;          % on line wavelength in nm
l_o2_off = 769.695;         % off line wavelength 
l_o2_abs_center = 769.795;  % absortption line center in nm
s0_o2 = .4889e-25;          % linestregth in cm/molecule
t0_o2 = 293.15;             % reference temperature in K
p0_o2 = 1;                  % reference pressure in atm
E_o2 = 1420.7631;           % ground state energy in 1/cm
gl0_o2 = 0.0312;            % halfwidth in 1/cm (gl0_o2*c = 1.085 GHz)
alpha_o2 = 0.63;

finesse_o2 = 2;
fsr_o2 = 0.1;  

% --- water vapor profile 
WV = wv_num';

%  ==================  Initial guess at temperature  ======================
T_g = T_surf - 6.5.*rkm;            % Assume a lapse rate of -6.5 K/km
T_g = T_g';                         % to stay consistant with DIAL_Temperature_Performance_v14_ksr.m program

rm = rm';                           % to stay consistant with DIAL_Temperature_Performance_v14_ksr.m program


%  using a linear fit to the temperature profile to get a surface
%  temperature and lapse rate.  This is useful in the second - fourth
%  iteration.
fitgamma1 = zeros(1,size(T_g,1));
fitintercept1 = zeros(1,size(T_g,1));
for i = 1:size(T_g,1)
    fitdata = fitlm(rm, T_g(i,:));
    fittable = fitdata.Coefficients;
    fitgamma1(i) = fittable{'x1', 'Estimate'};
    fitintercept1(i) = fittable{'(Intercept)', 'Estimate'};
end

%  =======================  Pressure profile  =============================
patm = P_surf'.*(fitintercept1'./(T_g)).^(.0341626./fitgamma1');
pkpa = 101.325.*patm;

% =========================================================================
%  Perturbative retrieval of alpha
% =========================================================================
%  -----------  zeroth order
oversample = movavg_r/37.5;             %averaged range size divided by the sample range bin size
i_range = size(rm, 2);
N_o2_on = o2on_avg(:,1:size(T,2));      % Number of time profiles may differ. Make T and N_o2 have the same.
N_o2_off = o2off_avg(:,1:size(T,2));

N_o2_on = N_o2_on';                     % to stay consistant with DIAL_Temperature_Performance_v14_ksr.m program
N_o2_off = N_o2_off';

ln_o2 = log(N_o2_on(:,1:end-oversample).*N_o2_off(:,1+oversample:end)./N_o2_on(:,1+oversample:end)./N_o2_off(:,1:end-oversample));
%delta_on = N_o2_on(:,1+oversample:end)./N_o2_on(:,1:end-oversample);
%delta_off = N_o2_off(:,1+oversample:end)./N_o2_off(:,1:end-oversample);

junk2 = N_o2_on(:,1:end-oversample).*N_o2_off(:,1+oversample:end)./N_o2_on(:,1+oversample:end)./N_o2_off(:,1:end-oversample);

%  ---  making the ln_o2 vectorplot(a the right light
i_ln_o2 = size(ln_o2, 2);
ln_o2_end = ln_o2(1, i_ln_o2);
ln_o2_end_v = repmat(ln_o2_end, 1, i_range - i_ln_o2);
ln_o2_total = [ln_o2  ln_o2_end_v];


alpha_0_ret = 1./2./rangebin.*ln_o2_total;
alpha_0_ret_filter = sgolayfilt(alpha_0_ret, 1, 1*oversample+1);
alpha_0_ret_smooth = smooth(alpha_0_ret_filter, 1*oversample+1)';
alpha_0_ret = alpha_0_ret_smooth;


%  -----------  first order
%  ---first calculate gi
%  spectral distributionusing the guessed at temperature profile
delfghz = -10:.06666666666666666:10;            %frequency scan in GHz
i_freq = size(delfghz, 2);
i_freq_middle = (i_freq+1)/2;
delfhz = delfghz.*1e9;
f_step = delfhz(1,2) - delfhz(1,1);
f_o2_on = c/(l_o2_on*1e-9);             %  center online frequency in hz
fscan_o2 = f_o2_on + delfhz;            %  frequency scan absolute
l_o2_scan = c./fscan_o2;                %  wavelength scan
wn_o2_scan = 1./l_o2_scan;              %  wavenumber scan in 1/m
wn_o2_scan_t = wn_o2_scan';
wno2_center = 1/(l_o2_on*1e-9);         %  center wavenumber in 1/m

middle = (size(wn_o2_scan, 2)+1)/2;     % gives me the index at the center of the frequency scan which is the operatiung wavelength

c_doppler_o2 = mo2*c^2/(8*pi*wno2_center^2*kb);



%  spectral distributionusing the guessed at temperature profile
doppler_o2_ret_un = ((c_doppler_o2./T_g).^.5).*exp(-c_doppler_o2.*(wno2_center-wn_o2_scan_t).^2./T_g);          
norm_o2_ret = trapz(doppler_o2_ret_un.*f_step);
doppler_o2_ret = doppler_o2_ret_un./norm_o2_ret; 

%  Checking to see if doppler_o2_ret is normalized to 1 when integrated across frequency
for i = 1:1:i_range
doppler_o2_ret_check(i) = trapz(doppler_o2_ret(:,i).*f_step);
end

%  ---  gi from retrieval data
%  ---  the modeled molecular backscatter
bmo2_ret = 374280*pkpa./T_g./l_o2_on.^4;
%  --- the aerosol backscatter coefficeint
ba_hsrl_ret_noise = betaa(:, pro_num)';
ba_hsrl_ret_filter = sgolayfilt(ba_hsrl_ret_noise, 3, 1*oversample+1);  %  Uses a Savitky-Golay filter%

beta_mn_o2_ret = bmo2_ret./(ba_hsrl_ret_filter+bmo2_ret);
beta_mn_o2_m_ret = repmat(beta_mn_o2_ret, i_freq,1);
gm_o2_ret = beta_mn_o2_m_ret.*doppler_o2_ret.*f_step;

beta_an_o2_ret = ba_hsrl_ret_filter./(ba_hsrl_ret_filter+bmo2_ret);
zero_m = zeros((i_freq_middle-1), size(rm, 2));
beta_an_o2_m_ret = [zero_m; beta_an_o2_ret; zero_m];
ga_o2_ret = beta_an_o2_m_ret;

g_o2_ret = ga_o2_ret+gm_o2_ret;

%  Checking to see if g_o2_ret is normalized to 1 when integrated across frequency
for i = 1:1:i_range
g_o2_ret_check(i) = trapz(g_o2_ret(:,i));
end

% ---  dg1/dr from retrieval data
dg_o2_ret = -1.*(g_o2_ret(:,1:size(rm, 2)-oversample)-g_o2_ret(:,1+oversample:size(rm, 2)));
drm = rm(1+oversample)-rm(1);
dg_dr_o2_ret = dg_o2_ret./drm;
dgdrlast = dg_dr_o2_ret(:,size(rm, 2)-oversample);
dgdrlast = repmat(dgdrlast, 1, oversample);

% ==========================================================================
dg_o2_ret2 = -1.*(g_o2_ret(:,1:size(rm, 2)-2*oversample)-g_o2_ret(:,1+2*oversample:size(rm, 2)));
drm2 = rm(1+2*oversample)-rm(1);
dg_dr_o2_ret2 = dg_o2_ret2./drm2;
dgdrlast2 = dg_dr_o2_ret2(:,size(rm, 2)-2*oversample);
dgdrlast2 = repmat(dgdrlast2, 1, 2*oversample);
dg_dr_o2_ret = [dg_dr_o2_ret2 dgdrlast2];
% =========================================================================

%dg_dr_o2_ret = [dg_dr_o2_ret dgdrlast];  %repeats the last column to make the matrix for dg/dr the same size as g

%  --- shifting by one rangebin
%dg_dr_o2_ret_unshift = dg_dr_o2_ret;
%dg_shift_short = dg_dr_o2_ret(:,oversample+1:size(rkm, 1));
%dg_end = dg_dr_o2_ret(:, size(rkm, 1));
%dg_end = repmat(dg_end, 1, oversample);
%dg_dr_o2_ret = [dg_shift_short dg_end];


%  ---  calculating the normalized molecular absorption lineshape absorption f
a_o2_ret = s0_o2.*(t0_o2./T_g);
c_o2_ret = exp(h.*c./kb.*(1./t0_o2-1./T_g).*E_o2);
sT_o2_ret = a_o2_ret.*c_o2_ret./100;                 % the divide by 100 makes the units m/molecule

gammaL_o2_ret = (gl0_o2*100).*(patm./p0_o2).*(t0_o2./T_g).^alpha_o2;        % times 100 makes the units 1/m
gammaD_o2_ret = (wno2_center/c).*(2.*kb.*log(2).*T_g./mo2).^.5;               % units are 1/m

j_wn_scan_o2 = size(wn_o2_scan,2);
shift_wn_o2 = 1/(l_o2_abs_center*1e-9) - 1/(l_o2_on*1e-9);
del_wn_o2 = wn_o2_scan - wno2_center - shift_wn_o2;

for i=1:1:size(rm, 2)      
    y = 0.8325546*gammaL_o2_ret(i)/gammaD_o2_ret(i);
    for j=1:1:j_wn_scan_o2;
        x = del_wn_o2(j)*.8325546/gammaD_o2_ret(i);
        dt = -10:.01:10;
        ft = exp(-dt.^2)./(y.*y+(x-dt).^2);
        convint = trapz(dt, ft);
        sigmaV_o2_ret(j,i) = sT_o2_ret(i)*0.12448*gammaL_o2_ret(i)/gammaD_o2_ret(i)^2*convint;      % units are m^2/molecule        
    end
end
linecenter_sigmaV = sigmaV_o2_ret((i_freq+1)/2,:);
linecenter_sigmaV_m = repmat(linecenter_sigmaV, i_freq, 1);
f_ret = sigmaV_o2_ret./linecenter_sigmaV_m;

% at line center, f_ret should have a value of 1.

%  --- Tmoth order
alpha_0_ret_m = repmat(alpha_0_ret, i_freq, 1);
alpha_0_freq_range = alpha_0_ret_m.*f_ret.*drm/oversample; 

for i = 1:1:size(rm, 2)
     for j = 1:1:i_freq
         integrand = alpha_0_freq_range(j,1:i);
         int = trapz(integrand);
         stuff(1,i) = int;
         Tmoth(j,i) = exp(-int);
     end
end

%  ---  O2 DIAL etalon
R_o2 = (((pi/finesse_o2/2)^2+1)^0.5-(pi/finesse_o2/2))^2;
fsr_o2_f = fsr_o2*1e-9*c/((l_o2_on*1e-9)^2);
theta_o2 = 2.*pi.*delfghz.*1e9./fsr_o2_f;
T_etalon_o2 = (1-R_o2).^2./(1+R_o2^2-2.*(R_o2).*cos(theta_o2));

T_etalon_o2 = T_etalon_o2';
T_etalon_o2_m = repmat(T_etalon_o2, 1, i_range);


eta = g_o2_ret.*T_etalon_o2.*Tmoth;
etals = eta.*(1-f_ret);
neta = dg_dr_o2_ret.*T_etalon_o2.*Tmoth;

eta2 = g_o2_ret.*T_etalon_o2;
neta2 = dg_dr_o2_ret.*T_etalon_o2;

for i=1:1:size(rm, 2)
    integrand1 = eta(:,i);
    integrand2 = etals(:,i);
    integrand3 = neta(:,i);
    integrand4 = eta2(:,i);
    integrand5 = neta2(:,i);
    eta_int(i) = trapz(integrand1);
    etals_int(i) = trapz(integrand2);
    neta_int(i) = trapz(integrand3);
    eta2_int(i) = trapz(integrand4);
    neta2_int(i) =trapz(integrand5);
end

%  --- first order perturbation
W1st = etals_int./eta_int;
G1st = neta_int./eta_int - neta2_int./eta2_int;
G1st_unshift = G1st;

G1_shift_short = G1st(1,.25*oversample+1:size(rkm, 1));
G1_end(1,1:.25*oversample) = G1st(1,size(rkm, 1));
G1st = [G1_shift_short G1_end];


alpha_1_ret = 0.5.*(alpha_0_ret.*W1st + G1st);

alpha_1_ret_filter = sgolayfilt(alpha_1_ret, 1, 1*oversample+1);
alpha_1_ret_smooth = smooth(alpha_1_ret_filter, 1*oversample+1)';
alpha_1_ret = alpha_1_ret_smooth;

alpha_0_1_ret = alpha_0_ret + alpha_1_ret;

alpha_1_ret_m = repmat(alpha_1_ret, i_freq,1);
alpha_1_freq_range = alpha_1_ret_m.*f_ret; 

%  ---  Tm1st

for i = 1:1:size(rm, 2)
     for j = 1:1:i_freq
         integrand = alpha_1_freq_range(j,1:i);
         int = drm*trapz(integrand);
         Tm1st(j,i) = exp(-int);
     end
end


%  ---  2nd order 
etat1 = eta.*(1-Tm1st);
netat1 = neta.*(1-Tm1st);
etat2 = eta.*(1-0);
netat2 = neta.*(1-0);
etat1f = eta.*(1-f_ret).*(1-Tm1st);

for i=1:1:size(rm, 2)
    integrand7 = etat1(:,i);
    integrand8 = netat1(:,i);
    integrand9 = etat1f(:,i);
    integrand10 = etat2(:,i);
    integrand11 = netat2(:,i);
    etat1_int(i) = trapz(integrand7);
    netat1_int(i) = trapz(integrand8);
    etat1f_int(i) = trapz(integrand9);
    etat2_int(i) = trapz(integrand10);
    netat2_int(i) = trapz(integrand11);
end

G2nd = (neta_int.*etat1_int./eta_int./eta_int - netat1_int./eta_int) - (neta_int.*etat2_int./eta_int./eta_int - netat2_int./eta_int);

G2nd_unshift = G2nd;

G2_shift_short = G2nd(1,.25*oversample+1:size(rkm, 1));
G2_end(1,1:.25*oversample) = G2nd(1,size(rkm, 1));%
G2nd = [G2_shift_short G2_end];


W2nd = etals_int.*etat1f_int./eta_int./eta_int - etat1f_int./eta_int;

alpha_2_ret = 0.5.*(alpha_1_ret.*W1st + alpha_0_ret.*W2nd + G2nd);

alpha_2_ret_filter = sgolayfilt(alpha_2_ret, 1, 1*oversample++1);
alpha_2_ret_smooth = smooth(alpha_2_ret_filter, 1*oversample+1)';
alpha_2_ret = alpha_2_ret_smooth;

alpha_o2_ret = alpha_0_ret + alpha_1_ret + alpha_2_ret;
alpha_o2_total_filter = sgolayfilt(alpha_o2_ret, 1, 1*oversample++1);
alpha_o2_ret = alpha_o2_total_filter;

%  ========================================================================
%  ================  Temperature Retrieval  ===============================
%  ========================================================================
loop = 20;                  % number of interations ror the temerpature retrieval 

T_ret(1,:) = T_g;
tg = T_g;

for lp = 2:1:loop+1
    % using a linear fit to the temperature profile to get a surface
    %  temeprature and lapse rate.  THis is useful in the second - fourth
    %  iteration.
    tgfit = tg(1,1:size(rm,2)-90);  %removes the furthest range bins befiorw fitting 
    rmfit = rm(1,1:size(rm,2)-90);
    fitdata = fitlm(rmfit, tgfit);
    fittable = fitdata.Coefficients;
    fitgamma1 = fittable{'x1', 'Estimate'};
    fitintercept1 = fittable{'(Intercept)', 'Estimate'};
    
    lapserate(lp-1) = fitgamma1;

    %  =======================  Pressure profile  =============================
    patm_g = psatm*(fitintercept1./(tg)).^(.0341626./fitgamma1);
    pg = patm_g(1,1)*101.325*1000;        %  Surface pressure in pascals
    pkpa_g = 101.325*patm_g;            %  Pressure profile in kPa 
    
    %  Find the new q term 
    Nair = No.*273.15.*patm_g./tg;    
    q = .2095.*(1-WV./Nair);

    % --- the lineshape value
    gammaL_o2_ret = (gl0_o2*100).*(patm_g./p0_o2).*(t0_o2./tg).^alpha_o2;        % times 100 makes the units 1/m
    gammaD_o2_ret = (wno2_center/c).*(2.*kb.*log(2).*tg./mo2).^.5;               % units are 1/m
    
    %  updating the lineshape function
    for i=1:1:size(rm, 2)       
       y = 0.8325546*gammaL_o2_ret(i)/gammaD_o2_ret(i);
       x = del_wn_o2(middle)*.8325546/gammaD_o2_ret(i);
       dt = -10:.01:10;
       ft = exp(-dt.^2)./(y.*y+(x-dt).^2);
       convint = trapz(dt, ft);
       f_o2_tg(i) = 0.12448*gammaL_o2_ret(i)/gammaD_o2_ret(i)^2*convint;      % units are m^2/molecule
    end
    
    %  calculating Delta T
    ep = E_o2*100*h*c;          % times 100 puts the units of E_o2 from 1/cm to 1/m
    exp1 = .03416262/fitgamma1;
    c1 = s0_o2/100*t0_o2*exp(ep/kb/t0_o2)*pg/kb/tsk^-exp1;
    c2 = tg.^(-exp1-2).*exp(-ep./kb./tg);
    c3 = (-exp1-2)./tg + ep./kb./tg./tg;

    deltaT = (alpha_o2_ret - c1.*q.*f_o2_tg.*c2)./(c1.*q.*f_o2_tg.*c2.*c3);
    deltaT_filter = sgolayfilt(deltaT, 1, 1*oversample+1);
    deltaT_smooth = smooth(deltaT_filter, 1*oversample+1)';

    
    %  ------- limit deltaT to plus or minus 5K  --------------------------
    for j = 1:1:size(rm, 2);
        if deltaT_smooth(1,j) > 5
            deltaT_smooth(1,j) = 5;
        end
        if deltaT_smooth(1,j) < -5
            deltaT_smooth(1,j) = -5;
        end
    end
    
    changeT(lp-1,:) = deltaT_smooth;
    new_T = tg+deltaT_smooth;
    T_ret(lp,:) = new_T;
    tg = new_T;
end 

T_ret_final = T_ret(loop, :);
T_ret_final_filter = sgolayfilt(T_ret_final, 1, 1*oversample+1);
T_ret_final_smooth = smooth(T_ret_final_filter, 1*oversample+1)';
T_ret_final = T_ret_final_smooth;