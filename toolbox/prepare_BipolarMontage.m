function prepare_BipolarMontage(Save,FileIn, FileOut)

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
% Authors: ?

[Root,file,ext]=fileparts(FileIn);
%Longitudinal bipolar montage
clear S
S.Filename=FileIn;
S.FileOut=FileOut;

if strcmp(Save,'False')
    Save=false;
else if strcmp(Save,'True')
        Save=true;
    end
end

S.Save=Save;
D = ImaGIN_BipolarMontage(S);

S.SaveFile=deblank(file);
disp(S.SaveFile(1:end-length(ext)-1));

end