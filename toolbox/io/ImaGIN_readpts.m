function [Name,Pos]=ImaGIN_readpts(filename)

fid=fopen(filename);
if fid~=-1
    tmp=fgetl(fid);
    tmp=fgetl(fid);
    Nelec=str2num(fgetl(fid));
    Name=cell(Nelec,1);
    Pos=zeros(Nelec,3);
    for i1=1:Nelec
        tmp=fgetl(fid);
        [Name{i1},tmp]=strtok(tmp);
        [tmp1,tmp]=strtok(tmp);
        Pos(i1,1)=str2num(tmp1);
        [tmp1,tmp]=strtok(tmp);
        Pos(i1,2)=str2num(tmp1);
        [tmp1,tmp]=strtok(tmp);
        Pos(i1,3)=str2num(tmp1);
    end
    fclose(fid);
else
    error('pts file not found')
end