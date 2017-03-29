function D = ImaGIN_Electrode(S)
% Set electrode positions

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
% Copyright (c) 2000-2017 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David

try
    t=S.Fname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

P=spm_str_manip(deblank(t(1,:)),'h');

try
    Position = S.Position;
catch
    try 
        filename = S.filenamePos;
        Position=load(filename);
    catch
        filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
        Position=load(filename);
    end
end


try
    Name = S.Name;
catch
    try
        filename = S.filenameName;
        Name = readName(filename,Position);
    catch
        filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
        Name = readName(filename,Position);
    end
end

try
    FileOut=S.FileOut;
catch
    S.FileOut = spm_input('Name of output file', '+1', 's');
end

try
    FileTxtOut = S.FileTxtOut;
catch
    disp('do not save the recording sensors')
end



for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    
    %     Cnames=chanlabels(D);
    Sensors=sensors(D,'EEG');
    Sensors.label=deblank(Sensors.label);
    
    if length(unique(Sensors.label))~=length(Sensors.label)
        warning('repeated label in the TRC file')
    end
    
    chan_not_assigned = cell(1,1);
    iter_bad = 1;
    for i1=1:length(Sensors.chantype)
        corresp=0;
        
        if strcmpi(chantype(D,i1),'eeg')
            for i2=1:length(Name{1})
                if strcmpi(Sensors.label{i1},Name{1}{i2})
                    Sensors.elecpos(i1,:)=Position(i2,:);
                    Sensors.chanpos(i1,:)=Position(i2,:);
                    corresp=corresp+1;
                end
                
                if strcmpi(Sensors.label{i1},Name{2}{i2})
                    Sensors.elecpos(i1,:)=Position(i2,:);
                    Sensors.chanpos(i1,:)=Position(i2,:);
                    corresp=corresp+1;
                end
            end
            
            Sensors.chantype{i1}='eeg';
            
            if corresp==0
                %try replacing ' by p
                for i2=1:length(Name{1})
                    if strcmpi(strrep(Sensors.label{i1},'''','p'),strrep(Name{1}{i2},'''','p'))
                        Sensors.elecpos(i1,:)=Position(i2,:);
                        Sensors.chanpos(i1,:)=Position(i2,:);
                        corresp=corresp+1;
                    end
                
                    if strcmpi(strrep(Sensors.label{i1},'''','p'),strrep(Name{2}{i2},'''','p'))
                        Sensors.elecpos(i1,:)=Position(i2,:);
                        Sensors.chanpos(i1,:)=Position(i2,:);
                        corresp=corresp+1;
                    end
                end
                
                
                if corresp == 0
                    warning([Sensors.label{i1} ' not assigned'])  %warning without identifier are never activated. won't see this message on most computers...
                    chan_not_assigned{iter_bad} = Sensors.label{i1};
                    iter_bad = iter_bad +1 ;
                end
%             elseif corresp>1
%                 warning(['Several electrodes corresponding to ' Sensors.label{i1}])
%                 disp(sprintf('Several electrodes corresponding to : %s', Sensors.label{i1}));
            end
            
        elseif strcmpi(chantype(D,i1),'ecg')
            Sensors.chantype{i1}='ecg';
        end
        
    end
    
   
   try
       fid = fopen(FileTxtOut,'w');
       recording_output = (D.sensors('eeg').label);
       fprintf(fid,'%s\n',recording_output{:});
       fclose(fid);
   catch
       disp('recording sensors not saved.')
   end
   
    D=sensors(D,'EEG',Sensors);
    newname=FileOut;
    D2=clone(D,newname, [D.nchannels D.nsamples D.ntrials]);
    D2(:,:,:)=D(:,:,:);
    save(D2);
    
end

end

    function Name = readName(filename,Position)
    fid=fopen(filename);
    for i1=1:size(Position,1)
        Name{i1}=fgetl(fid);
    end
    fclose(fid);
    end
