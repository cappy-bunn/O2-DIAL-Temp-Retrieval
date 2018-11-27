oid = 'raob_soundings24581.cdf';

f = ncinfo(oid);
nvars = length(f.Variables);
for k = 1:nvars
   varname=f.Variables(k).Name;
   disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' oid ''',''' varname ''');']);
end

% Time/Date
synTime = datetime(synTime,'convertfrom','posixtime');  % Convert from unix time to normal time/date

% === Set mask to NaNs ===

% Mandatory Levels
tpMan(tpMan == max(max(tpMan))) = NaN;
tpMan(tpMan == max(max(tpMan))) = NaN;
prMan(prMan == max(max(prMan))) = NaN;
htMan(htMan == max(max(htMan))) = NaN;
htMan(1,:) = 0;     % Set station location to zero instead of elevation

% Significant Levels
tpSigT(tpSigT == max(max(tpSigT))) = NaN;
htSigT(htSigT == max(max(htSigT))) = NaN;


% === Interleave mandatory and significant observations to get full profiles ===

% Select time and date
time_date = '29-Oct-2018 00:00:00';
[~,t_profile] = min(abs(synTime-time_date));

% (specific ONLY for Oct 29, 2018 at 0000 UTC)
temp_oct29_00UTC = [tpMan(1,t_profile);tpMan(5,t_profile);tpSigT(1,t_profile);tpMan(6,t_profile)];
height_oct29_00UTC = [htMan(1,t_profile);htMan(5,t_profile);htSigT(1,t_profile);htMan(6,t_profile)];
figure
plot(temp_oct29_00UTC,height_oct29_00UTC)

% Select time and date
time_date = '29-Oct-2018 12:00:00';
[~,t_profile] = min(abs(synTime-time_date));

% (specific ONLY for Oct 29, 2018 at 1200 UTC)
temp_oct29_12UTC = [tpMan(1,t_profile);tpSigT(1:3,t_profile);tpMan(5,t_profile);tpSigT(4:5,t_profile);tpMan(6,t_profile);tpSigT(6,t_profile)];
height_oct29_12UTC = [htMan(1,t_profile);htSigT(1:3,t_profile);htMan(5,t_profile);htSigT(4:5,t_profile);htMan(6,t_profile);htSigT(6,t_profile)];
figure
plot(temp_oct29_12UTC,height_oct29_12UTC)

