oid = 'air.2018.nc';

f = ncinfo(oid);
nvars = length(f.Variables);
for k = 1:nvars
   varname=f.Variables(k).Name;
   disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' oid ''',''' varname ''');']);
end

oid2 = 'hgt.2018.nc';

f2 = ncinfo(oid2);
nvars2 = length(f2.Variables);
for k = 1:nvars2
   varname2=f2.Variables(k).Name;
   disp(['Reading:  ' varname2]);
   eval([varname2 ' = ncread(''' oid2 ''',''' varname2 ''');']);
end


% Select time and date (measurements 4 times daily: 0000,0600,1200,1800)
t_profile = '29-Oct-2018 00:00:00';

% Select location latitude and longitude
lat_NCAR = 40.0378;                     % Degrees N
lon_NCAR = 360 - 105.2416;              % Degrees E
[~,lat_ind] = min(abs(lat-lat_NCAR));   % Find index of NCAR lat in lat matrix
[~,lon_ind] = min(abs(lon-lon_NCAR));   % Find index of NCAR lon in lon matrix

air_NCAR = squeeze(air(lon_ind,lat_ind,:,:));   % Air temp profile at NCAR in 2018
hgt_NCAR = squeeze(hgt(lon_ind,lat_ind,:,:));   % Geopotential height at NCAR in 2018

% Time is hours since 01-Jan-1800, 0000UTC
t1 = datetime('01-Jan-2018 00:00:00');
t2 = datetime(t_profile);
[h,m,s] = hms(t2 - t1);

t = time(1) + h;                    % Convert to hrs since 01-Jan-1800 0000UTC
[hr_diff,hr_ind] = min(abs(time-t));

air_NCAR_sp = air_NCAR(:,hr_ind);   % Single air temp profile at time t
hgt_NCAR_sp = hgt_NCAR(:,hr_ind);   % Single geopotential height profile at time t

max_hgt = 6000;                     % Max height to plot [m]
[~,max_hgt_ind] = min(abs(hgt_NCAR_sp-max_hgt));

% Height vs Air Temp plot at NCAR at time t
plot(air_NCAR_sp(1:max_hgt_ind),hgt_NCAR_sp(1:max_hgt_ind))
title(sprintf('NCEP Air Temperature at NCAR on %s UTC', t_profile));
xlabel('Temperature (K)')
ylabel('Range (m)')

% Super rough estimate of lapse rate
lapserate = (air_NCAR_sp(max_hgt_ind) - air_NCAR_sp(1))/((hgt_NCAR_sp(max_hgt_ind) - hgt_NCAR_sp(1))/1000);
