function evt=ImaGIN_ReadXltekEvents(filename)


%Read Xltek events (.txt file separated from edf file with the same name)

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
