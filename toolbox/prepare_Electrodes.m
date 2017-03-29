function prepare_Electrodes(source, polarity, FileIn, implantationFile, FileOut, FileTxtOut)
% Add electrode position
%
%source : 'classic' (a file with the electrodes names and another one with
%the MNI positions), implantationFile is the source directory of the 2
%files
%      or 'intranat' (a file -implantationFile- with electrodes names, position, atlas,...)

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
%
% Authors: ?

[Root,file,~] = fileparts(FileIn);

clear S
S.Fname = fullfile(Root,[spm_str_manip(file,'r') '.mat']);

if strcmp(source,'classic')
    S.filenamePos = fullfile(implantationFile,'Electrodes_Pos_MNI.txt');
    S.filenameName = fullfile(implantationFile,'Electrodes_Name.txt');
end

%for files generated with Intranat
if strcmp(source, 'intranat')
    data=extractIntranatFile(implantationFile);
    
    %identification of the line with the MNI positions
    for i1=1:size(data,2)
        for i2=1:size(data{i1},1)
            if strcmp(data{i1}{i2},'MNI')
                indi=i1;
                indj=i2;
            end
        end
    end
    positiontmp = data{indi};
    nametmp = data{1};
    Position = [];
    Name = cell(2,1);
    if strcmp(polarity,'monopolar')
        monop=[];
        for i=(indj+1):length(nametmp)
            if isempty(strfind(nametmp{i},'-'))
                ind1 = regexp(lower(nametmp{i}),'[a-z\'']');
                ind2 = regexp(lower(nametmp{i}),'[0-9]');
                num = num2str(str2num(nametmp{i}(ind2)));
                num2 = nametmp{i}(ind2);
                name = [lower(nametmp{i}(ind1)) num];
                name2 = [lower(nametmp{i}(ind1)) num2];
                Name{1} = cat(1,Name{1},{name});
                Name{2} = cat(1,Name{2},{name2});
                Position = cat(1,Position,str2num(positiontmp{i}));
            end
        end
    else
        if strcmp(polarity, 'bipolar')
        bip=[];
        for i=4:length(nametmp)
            if ~isempty(strfind(nametmp{i},'-')) 
                ind1 = regexp(lower(nametmp{i}),'[a-z\'']');
                sep = regexp(nametmp{i},'-');
                ind2 = regexp(nametmp{i},'\d');
                if ~isempty(ind1) && ~isempty(ind2)
                    name = [lower(nametmp{i}(ind1 < sep)) num2str(str2num(nametmp{i}(ind2(ind2 < sep)))) lower(nametmp{i}(ind1 < sep)) num2str(str2num(nametmp{i}(ind2(ind2 > sep))))];
                    name2 = [lower(nametmp{i}(ind1 < sep)) nametmp{i}(ind2(ind2 < sep)) lower(nametmp{i}(ind1 < sep)) nametmp{i}(ind2(ind2 > sep))];
                    Name{1} = cat(1,Name{1},{name});
                    Name{2} = cat(1,Name{2},{name2});
                    Position = cat(1,Position,str2num(positiontmp{i}));
                end
            end
        end
        end
    end
S.Position = Position;
S.Name = Name;
end


S.FileTxtOut = FileTxtOut;
S.FileOut = FileOut;
D = ImaGIN_Electrode(S);

end


function data = extractIntranatFile(filename)
delimiter = '\t';
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');
data = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
end