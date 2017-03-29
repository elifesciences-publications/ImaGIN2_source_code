function D = ImaGIN_causality_average(S)
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

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TF',0);
try
    D1 = S.D1;
catch
    D1 = spm_select(inf, '\.mat$', 'Select causality data file',[],pwd,'ca1_');
end

P = spm_str_manip(D1, 'H');

try
    if size(D1,1)==1
        D1 = spm_eeg_load(D1);
        d=spm_str_manip(D1,'h');
        e=spm_str_manip(D1,'t');
        tmp=strfind(e,'_');
        e=[e(1:tmp(1)+2) '2' e(tmp(1)+4:end)];
        t=fullfile(d,e);
        D2=spm_eeg_ldata(t);    %time lag
    else
        DD=D1;
        clear D1
        for i1=1:size(DD,1)
            D1{i1}= spm_eeg_load(deblank(DD(i1,:)));
            d=spm_str_manip(deblank(DD(i1,:)),'h');
            e=spm_str_manip(deblank(DD(i1,:)),'t');
            tmp=strfind(e,'_');
            e=[e(1:tmp(1)+2) '2' e(tmp(1)+4:end)];
            t=fullfile(d,e);
            D2{i1}=spm_eeg_load(t);    %time lag
        end
    end
catch
    error(sprintf('Trouble reading file %s', D));
end

if iscell(D1)
    DD1=D1;clear D1
    DD2=D2;clear D2
    data1=0;
    data2=0;
    for i1=size(DD1,1)
        data1=data1+double(DD1{i1}(:,:,:));
        data2=data2+double(DD2{i1}(:,:,:));
    end
    data1=data1./size(DD1,1);
    data2=data2./size(DD2,1);
    Name=spm_input('Name of new file', '+1', 's');
    D1=DD1{1};
    D1=rmfield(D1,'data');
    D2=DD2{1};
    D2=rmfield(D2,'data');
    if isempty(Name)%Assume namefiles are numbered, have the same events
        for i1=1:length(D1.fnamedat)
            if ~strcmp(DD1{1}.fnamedat(i1),DD1{2}.fnamedat(i1))
                i2=find(D1.fnamedat(i1+1:end)=='_');
                D1.fnamedat=[D1.fnamedat(1:i1-1) 'Mean_' D1.fnamedat(i1+i2(1)+1:end)];
                break
            end
        end
        for i1=1:length(D2.fnamedat)
            if ~strcmp(DD2{1}.fnamedat(i1),DD2{2}.fnamedat(i1))
                i2=find(D2.fnamedat(i1+1:end)=='_');
                D2.fnamedat=[D2.fnamedat(1:i1-1) 'Mean_' D2.fnamedat(i1+i2(1)+1:end)];
                break
            end
        end
    else
        D1.fnamedat=[Name '_ca1.dat'];
        D2.fnamedat=[Name '_ca2.dat'];
    end
    D1.fname=[D1.fnamedat(1:end-3) 'mat'];
    D2.fname=[D2.fnamedat(1:end-3) 'mat'];
    P=deblank(P(1,:));
    fpd = fopen(fullfile(P, D1.fnamedat), 'w');
    for i=1:D1.Nevents;
        D1.scale(:, i) = spm_eeg_write(fpd, squeeze(data1(:,:,i)), 2, D1.datatype);
    end
    fclose(fpd);
    fpd = fopen(fullfile(P, D2.fnamedat), 'w');
    for i=1:D2.Nevents;
        D2.scale(:, i) = spm_eeg_write(fpd, squeeze(data2(:,:,i)), 2, D2.datatype);
    end
    fclose(fpd);
    cd(D1.path)
    if str2num(version('-release'))>=14
        D=D1;
        save(fullfile(P, D1.fname), '-V6', 'D');
        D=D2;
        save(fullfile(P, D2.fname), '-V6', 'D');
    else
        D=D1;
        save(fullfile(P, D1.fname), 'D');
        D=D2;
        save(fullfile(P, D2.fname), 'D');
    end
else
    error('Select more than one file');
end
