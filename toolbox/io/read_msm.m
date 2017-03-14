function [data header] = read_msm(filename)
% Import msm and msr files
% Returns data  = eeg data where each column is a channel
%         header = information about channel names, sampling rate, etc
%        
% Author: Cristian Donos 
% Physics Dept, University of Bucharest, March 2015
% email: cristian.donos@g.unibuc.ro
%


%% header - read information about the contacts
fid = fopen([filename '.msr']);

str = fgetl(fid);
header.record_length = str2num(str(strfind(str,'=')+1:end));

str = fgetl(fid);
header.numchannels = str2num(str(strfind(str,'=')+1:end));

str = fgetl(fid);
header.unitmeas =  strtrim(str(max(find(isspace(str))):end));

str = fgetl(fid);
header.unittime = strtrim(str(max(find(isspace(str))):end));

str = fgetl(fid);
str = str(max(find(isspace(str))):end);
if strcmp(header.unittime,'ms')
     header.samplingrate(1:header.numchannels) = 1000 / str2num(str(strfind(str,'(')+1:end-1)) - str2num(str(1:strfind(str,'(')-1));
else
     disp(' MSM file import not available for other units than "ms"');
     return
end

% contact labels
str = fgetl(fid); % Labels lines start with the next fgetl
header.channels = [];
ix = 1;
str = fgetl(fid);
[tok remaining ] = strtok(str);
header.channels{ix} = tok;
ix = ix+1;
while ~isempty(remaining)
    [tok remaining ] = strtok(remaining);
    header.channels{ix} = tok;
    ix = ix+1;      
end
header.channels(end)=[]; % delete last cell as it has been created by the last ix = ix+1
header.channels = header.channels';

fclose(fid);    

%% data - read eeg data
fid = fopen([filename '.msm']);
data = fread(fid,Inf,'float32');
data =  reshape(data,[header.numchannels length(data)/header.numchannels ])';
fclose(fid);

end
