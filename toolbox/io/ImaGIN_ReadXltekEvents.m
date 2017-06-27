function evt=ImaGIN_ReadXltekEvents(filename)
% Read Xltek events (.txt file separated from edf file with the same name)

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
% Authors: Olivier David, 2017

fid=fopen(filename);
if fid~=-1
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    
    %read start time
    [tmp1,tmp2]=strtok(tmp);
    starttime=strtok(tmp2);
    starttime=convertHHMMSStoS(num2str(starttime));
    
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    
    evt=[];
    N=0;
    ok=1;
    while ok
        tmp=fgetl(fid);
        if tmp==-1
            ok=0;
        else
            [tmp1,tmp2]=strtok(tmp);
            [tmp1,tmp2]=strtok(tmp2);
            time=convertHHMMSStoS(num2str(tmp1));
            comment=strtrim(tmp2);
            switch comment
                case {'XLEvent','XLSpike','Gain/Changement de filtre','e','h','s','Opened relay'}
                otherwise
                    if length(comment)>12
                        if strcmp(comment(1:12),'Closed relay')
                        else
                            N=N+1;
                            evt(N).type  = comment;
                            evt(N).time = time-starttime;
                            evt(N).value= N;
                        end
                    else
                        N=N+1;
                        evt(N).type  = comment;
                        evt(N).time = time-starttime;
                        evt(N).value= N;
                    end
            end
                        
        end
    end
    fclose(fid);
else
    error('event file not found')
end

end


%convert a time in HHMMSS format in seconds
function H=convertHHMMSStoS(hms)
    h=str2num(hms(1:2));
    m=str2num(hms(4:5));
    s=str2num(hms(7:8));

    H=h*3600+m*60+s;
end

