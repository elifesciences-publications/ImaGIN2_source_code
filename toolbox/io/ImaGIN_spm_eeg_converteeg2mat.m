function D = ImaGIN_spm_eeg_converteeg2mat(S)
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
% $Id: spm_eeg_converteeg2mat.m 193 2005-07-01 13:37:59Z james $


try
    S2.dataset=S.dataset;
catch
    [S.dataset, sts] = spm_select(1, '.*', 'Select M/EEG data file');
    S2.dataset=S.dataset;
    if ~sts, return; end
end
Filename=S2.dataset(1:end-4);

fmt=lower(spm_str_manip(S.dataset,'e'));


% which species?
try
    S2.Atlas = S.Atlas;
catch
    %     Ctype = {
    %             'Human',...
    %             'Rat',...
    %             'Mouse'};
    %     str   = 'Select atlas';
    % 	Sel   = spm_input(str, '+1', 'm', Ctype);
    % 	S2.Atlas = Ctype{Sel};
    S2.Atlas = spm_input('Select atlas', '+1','Human|Rat|Mouse');
end


try
    S2.FileOut=S.FileOut;
catch
    S2.FileOut = spm_str_manip(S.dataset(1:end-4),'t');
end

% % which channel template file
% try
%     Fchannels = S.Fchannels;
% catch
%     Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
% end
% S2.Fchannels = Fchannels;

Nfiles = size(S.dataset, 1);
Cnames ='';
switch fmt
    case {'smr'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Fchannels;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.bipole = S.Montage;
            end
            try
                S2.epochlength = S.epochlength;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.channel = S.channel;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_spike2_mono(S2);
        end

    case {'eeg'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Fchannels;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.bipole = S.Montage;
            end
            try
                S2.epochlength = S.epochlength;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.channel = S.channel;
            end
            try
                S2.SaveFile = S.SaveFile;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_elan(S2);
        end

    case {'asc','txt'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Fchannels;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.channel = S.channel;
            end
            try
                S2.Radc = S.Radc;
            end
            try
                S2.Nevent = S.Nevent;
            end
            try
                S2.filenamePos=S.filenamePos;
            end
            try
                S2.filenameName=S.filenameName;
            end
            try
                S2.MontageName = S.MontageName;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_ascii(S2);
        end

    case {'trc'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.pts = S.pts;
            end
            try
                S2.System = S.System;
            end
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Fchannels;
            end
            try
                S2.channel = S.channel;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.Montage = S.Montage;
            end
            try
                S2.MontageName = S.MontageName;
            end
            try
                S2.NeventType = S.NeventType;
            end
            try
                S2.event =  S.event;
            end
            try
                S2.event_file=S.event_file;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.bipole=S.bipole;
            end
            try
                S2.loadevents=S.loadevents;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_micromed_mono(S2);
        end

    case {'msm'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.Atlas = S.Atlas;
            end
            try
                S2.SEEG = S.SEEG;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.SaveFile=S.SaveFile;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_msm_mono(S2);
        end

    case {'bin'}
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.bipole = S.Bipole;
            end
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Montage;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.Nevent=S.Nevent;
            end
            try
                S2.SEEG=S.SEEG;
            end
            try
                S2.filenamePos=S.filenamePos;
            end
            try
                S2.filenameName=S.filenameName;
            end
            try
                S2.MontageName = S.MontageName;
            end
            try
                S2.SaveFile=S.SaveFile;
            end
            D{i1} = ImaGIN_spm_eeg_rdata_deltamedbin_mono(S2);
        end

       
    case {'edf'}
        
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
            try
                S2.channel = S.channel;
            end
            try
                S2.Atlas = S.Atlas;
            end
            try
                S2.SEEG = S.SEEG;
            end
            try
                S2.Bipolar=S.Bipolar;
            end
            try
                S2.coarse = S.coarse;
            end
            try
                S2.SizeMax = S.SizeMax;
            end
%             try
%                 S2.MontageName=S.MontageName;
%             end
%             try
%                 S2.filenamePos=S.filenamePos;
%             end
%             try
%                 S2.Fchannels=S.Fchannels;
%             end
%             try
%                 S2.CreateTemplate = S.CreateTemplate;
%             end

            spm('Pointer','Watch');
            D = ImaGIN_spm_eeg_rdata_edf(S2);
        end
        spm('Pointer','Arrow');
        
    case {'e'}
        
        for i1 = 1:Nfiles
            S2.Fdata = deblank(S.dataset(i1, :));
%             try
%                 S2.Bipolar=S.Bipolar;
%             end
%             try
%                 S2.bipole = S.Bipole;
%             end
            try
                S2.CreateTemplate = S.CreateTemplate;
            end
            try
                S2.Fchannels = S.Montage;
            end
            try
                S2.coarse = S.coarse;
            end
%             try
%                 S2.Nevent=S.Nevent;
%             end
            try
                S2.SEEG=S.SEEG;
            end
            try
                S2.filenamePos=S.filenamePos;
            end
            try
                S2.filenameName=S.filenameName;
            end
            try
                S2.MontageName = S.MontageName;
            end
            try
                S2.SaveFile=S.SaveFile;
            end
            try
                S2.channel = S.channel;
            end
            D = ImaGIN_spm_eeg_rdata_nicolet_mono(S2);
        end
        spm('Pointer','Arrow');
        
    otherwise
        error('Unknown format');
end

end


