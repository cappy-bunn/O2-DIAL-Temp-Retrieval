oid = 'hgt.2018.nc';

f = ncinfo(oid);
nvars = length(f.Variables);
for k = 1:nvars
   varname=f.Variables(k).Name;
   disp(['Reading:  ' varname]);
   eval([varname ' = ncread(''' oid ''',''' varname ''');']);
end