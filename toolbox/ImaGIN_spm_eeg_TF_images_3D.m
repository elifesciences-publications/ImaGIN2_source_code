function D = ImaGIN_spm_eeg_TF_images_3D(S)

%make analyse images for intracerebral EEG
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
	D = S.D;
catch
	D = spm_select(inf, '\.mat$', 'Select EEG mat file');
	
end
P = spm_str_manip(D, 'H');

try
	D = spm_eeg_load(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

if isfield(D, 'Nfrequencies') & ndims(D)==3
	try
		fmt = S.fmt;
	catch
		spm_input('average over ...', 1, 'd')
		Ctype = {
			'electrodes',...
				'frequency'};
		str   = 'Average over which dimension';
		Sel   = spm_input(str, 2, 'm', Ctype);
		fmt = Ctype{Sel};
	end
	
	switch fmt
		case {'electrodes'}
			try
				D.electrodes_of_interest = S.thresholds.elecs;
			catch 
				str = 'electrodes[s]';
				Ypos = -1;
				
				while 1
					if Ypos == -1   
						[D.electrodes_of_interest, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
					else
						D.electrodes_of_interest = spm_input(str, Ypos, 'r', [], [1 Inf]);
					end
					
					
					t=1:D.nchannels;
					tmp=[];
					for en=D.electrodes_of_interest;
						if isempty(find(t==en))
							tmp=[tmp,en];
						end
					end
					if isempty(tmp) break, end
				end
			end
			try
				D.Nregion = S.region_no;
			catch 
				str = 'region number';
				Ypos = -1;
				
				while 1
					if Ypos == -1   
						[D.Nregion, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
					else
						D.Nregion = spm_input(str, Ypos, 'r', [], [1 Inf]);
					end
					if ~isempty(D.Nregion) break, end
					str = 'No data';
				end
				
            end
            Types={};
            for i1=1:length(Events)
                trouve=0;
                for i2=1:length(Types)
                    if ~strcmp(Types{i2},Events(i1).type) && i2==length(Types) && ~trouve
                        Types{end+1}=Events(i1).type;
                    elseif strcmp(Types{i2},Events(i1).type)
                        trouve=1;
                    end
                end
                if isempty(Types)
                    Types{1}=Events(1).type;
                end
            end
            for i = 1 : length(Types)
                Itrials = find(strcmp(Events.type == Types(i)));
                cd(D.path)
                dname = sprintf('%dROI_TF_trialtype%d', D.Nregion, Types(i));
                [m, sta] = mkdir(dname);
                cd(dname);
				
				for l = Itrials
					% if single trial data make new directory for single trials,
					% otherwise just write images to trialtype directory
					if size(D,4) ~= length(Types)
						% single trial data
						dname = sprintf('trial%d.img', l);
						fname = dname;
						[m, sta] = mkdir(dname);
						cd(dname);
					else
						fname = 'average.img';
					end
					data=squeeze(mean(D(D.electrodes_of_interest,:,:,i),1));	
					V.fname = fname;
					V.dim = [D.Nfrequencies D.nsamples  1 ];
					V.dt=[spm_type('float64') 0]; %%%check later with john
					V.mat = eye(4);
					V.pinfo = [1 0 0]';
					
					spm_write_vol(V, data); % d is data
				end
				
			end
			
        case {'frequency'}
            try
                D.Frequency_window = S.freqs;
                Ypos = -1;
                while 1
                    if Ypos == -1
                        Ypos = '+1';
                    end
                    
                    inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2))
                    if ~isempty(inds) break, end
                    str = 'No data in range';
                end
            catch 
                str = 'Frequency window';
				
				Ypos = -1;
				while 1
					if Ypos == -1
						Ypos = '+1';
					end
                    [D.Frequency_window, Ypos] = spm_input(str, Ypos, 'r', [], 2);
                    
                    inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2))
                    if ~isempty(inds) break, end
                    str = 'No data in range';
                end
            end
            data=squeeze(mean(D(:,inds,:,:),2));
            D=clone(D,['F' num2str(D.Frequency_window(1)) '_' num2str(D.Frequency_window(2)) '_' D.fnamedat], [size(D,1) size(D,3) 1]);
            D(:,:)=data;
            
           % Ctf = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
           % D.channels.Bad = setdiff(1:Ctf.Nchannels,D.tf.channels);
            D.time=D.tf.time;
            D=rmfield(D,'Nfrequencies');
            save(D);
            S.Fname=D.fname;
            
            try
                n = S.n;
            catch
                S.n = spm_input('Output image spatial resolution [mm]', '+1', 'n', '3', 1);
                n=S.n;
            end
            
            if length(n) > 1
                error('n must be scalar');
            end
            
            try
                interpolate_bad = S.interpolate_bad;
            catch
                S.interpolate_bad = spm_input('Interpolate bad channels or mask out?',...
                    '+1', 'b', 'Interpolate|Mask out', [1,0]);
            end
            ImaGIN_spm_eeg_convertmat2ana_3D(S);
    end
else
    clear S;
    S.Fname = fullfile(D.path, D.fname);
    ImaGIN_spm_eeg_convertmat2ana_3D(S);
end
