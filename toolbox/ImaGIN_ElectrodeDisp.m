function ImaGIN_ElectrodeDisp(varargin)
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
% Authors: Francois Tadel, 2017

if (nargin==1)
    if varargin{1}==1
        [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Create montage',0);
        Action='Create';
    elseif varargin{1}==2
        [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display montage',0);
        Action='Display';
    end
elseif nargin==0
    Action = spm_input('Montage ',1,'Create|Display');
end

if strcmp(Action,'Create')
    ImaGIN_MontageCreate;
elseif strcmp(Action,'Display')
    ImaGIN_MontageDisplay;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_MontageCreate
    tmp = spm_select(1, '\.txt$', 'Select file for names of electrodes');
    if isempty(tmp)
        return;
    end
    fid=fopen(tmp);
    n=0;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        n=n+1;
        Cnames{n,1}=tline;
    end
    fclose(fid);

    tmp = spm_select(1, '\.txt$', 'Select file for positions of electrodes');
    Cpos2=load(tmp);
    if size(Cpos2,1)~=3
        Cpos2=Cpos2';
    end
    Cpos=Cpos2(1:2,:);
    Cpos=Cpos-min(Cpos(:));
    Cpos=Cpos./max(Cpos(:));

    if size(Cpos,2)==length(Cnames)
        Nchannels=size(Cpos,2);
    else
        error('Number of electrodes positions different from number of electrodes names')
    end

    Species = spm_input('Create template ','+1','Rat|Mouse|Human');
    switch Species
        case{'Rat','Mouse'}
            Cpos2=10*Cpos2;
    end

    Files=what(fullfile(spm('dir'),'EEGtemplates'));

    ok=1;
    while ok
        Filename=spm_input('Name of montage', '+1','s');
        ok=0;
        for i1=1:length(Files.mat)
            if strcmp(lower(Files.mat{i1}),lower(Filename))
                spm_input('Specified file already exists', '+1','d')
                ok=1;
                break
            end
        end
    end
    Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
    save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_MontageDisplay
    FS1 = spm('FontSize', 14);
    Action = spm_input('Montage ','+1','2D|3D');
    Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
    M=load(Fchannels);

    switch Action
        case{'2D'}
            switch M.Species
                case{'Rat'}
                    Im=load(fullfile(spm('dir'), 'EEGtemplates','BregmaZoom'));
            end
            Fgraph  = spm_figure('GetWin','Graphics');

            figure(Fgraph);clf

            imagesc(Im.x,Im.y,Im.Bregma)
            axis equal
            axis tight
            title('Paxinos coordinates in mm','FontSize',FS1)

            hold on
            for i1=1:M.Nchannels
                c=plot(M.Cpos2(1,i1),M.Cpos2(2,i1),'r.','Markersize',40);
                set(c,'Tag',sprintf('%d',i1))
                h=text(M.Cpos2(1,i1)+.5,M.Cpos2(2,i1),M.Cnames{i1});
                set(h,'FontWeight','Bold','FontSize',FS1,'Color','b','HorizontalAlignment','left');
            end
            hold off
            colormap('gray')

        case{'3D'}
            sdip.n_seeds = 1;
            sdip.n_dip = M.Nchannels;
            sdip.Mtb = 1;
            sdip.j{1} = zeros(3*M.Nchannels, 1);
            sdip.loc{1} = M.Cpos2;
            sdip.Names=M.Cnames;

            tmp=spm_input('T1 template overlay ','+1','Yes|No');
            if strcmp(tmp,'Yes')
                switch M.Species
                    case{'Human'}
                        P = fullfile(spm('dir'),'canonical','avg152T1.nii');
                    case{'Rat'}
                        P=fullfile(spm_str_manip(spm('dir'),'h'),'ratlas5','template','template_T1.img');
                    case{'Mouse'}
                        P=fullfile(spm_str_manip(spm('dir'),'h'),'mousatlas5','template','template_T1.img');
                end
            else
                P   = spm_select(1,'image','Select image for rendering on');
            end
            ImaGIN_DrawElectrodes('Init', sdip,P)
    end
end


