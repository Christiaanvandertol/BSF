function [wl, data, time] = readFXBox(filename)

fileID = fopen(filename);
T = readtable(filename, 'TreatAsEmpty', '#N/D');
fclose(fileID);
x = table2array(T);

format = '%4c';
for n = 1:size(x,2)
    format =[format '%2c%d%c%d%c%d%c']; %#ok<AGROW>
end

fid = fopen(filename);
line = fgetl(fid);
z = sscanf(line,format);
fclose('all');
hr = z(7:8:end);
minute = z(9:8:end);
sec = z(11:8:end);
wl = x(:,1);
%data = x(:,2:Last+1);
data = x(:,2:end);
time = datenum(0,0,0,hr,minute,sec); %#ok<*SAGROW>