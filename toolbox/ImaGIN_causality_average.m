function D = ImaGIN_causality_average(S)

%Average TF in a given frequency band (from 3D to 2D data)
%Olivier David

% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_converteeg2mat(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be either 'CNT' or 'BDF'
%    Mname     - char matrix of input file name(s)
%    Fchannels - String containing name of channel template file
%_______________________________________________________________________
%
% spm_eeg_converteeg2mat is a user interface to convert EEG-files from their
% native format to SPM's data format. This function assembles some
% necessary information before branching to the format-specific conversion
% routines.
% The user has to specify, by either using struct S or the GUI, a 'channel
% template file' that contains information about the (approximate) spatial
% positions of the channels.
%
% Output: The converted data are written to files. The header
% structs, but not the data, are returned in D as a cell vector of structs.
%_______________________________________________________________________
%
% Additional formats can be added by (i) extending the code below in a
% straightforward fashion, (ii) providing a new channel template file and
% (iii) adding the actual conversion routine to the SPM-code.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_TF_images.m 299 2005-11-15 15:25:17Z james $

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
