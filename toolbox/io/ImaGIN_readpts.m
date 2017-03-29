function [Name,Pos]=ImaGIN_readpts(filename)
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
% Authors: Olivier David

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