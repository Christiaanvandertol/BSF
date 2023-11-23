function [wl, data, time] = readFXBox(filename,numchar)

fileID = fopen(filename);
T = readtable(filename, 'TreatAsEmpty', '#N/D');
fclose(fileID);
x = table2array(T);

if nargin<2
    numchar = 4;
end

format = [ '%' num2str(numchar) 'c'];
for n = 1:size(x,2)
%    format =[format '%2c%d%c%d%c%d%c']; %#ok<AGROW>
    format =[format '%3c%d%c%d%c%f%c']; %#ok<AGROW>
end

fid = fopen(filename);
line = fgetl(fid);
z = sscanf(line,format);
fclose('all');
%hr = z(7:8:end);
%minute = z(9:8:end);
%sec = z(11:8:end);

hr = z(8-4+numchar:9:end);
minute = z(10-4+numchar:9:end);
sec = z(12-4+numchar:9:end);

wl = x(:,1);
%data = x(:,2:Last+1);
data = x(:,2:end);

time = datenum(0,0,0,hr,minute,sec); %#ok<*SAGROW>

if abs(length(time)- size(data,2))>0
    data = data(:,1:length(time));        
end
