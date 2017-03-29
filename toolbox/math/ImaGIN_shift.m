function J	= ImaGIN_shift(J,pas,dim)
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

j	= J;
if min(size(J)) == 1
   if pas > 0
      J(1:pas)	= j(length(J)-pas+1:length(J));
      J(pas+1:length(J))	= j(1:length(J)-pas);
   elseif pas < 0
      J(1:length(J)+pas)	= j(-pas+1:length(J));
      J(length(J)+pas+1:length(J))	= j(1:-pas);
   end
elseif nargin == 3
    if ndims(J) == 2
        if dim == 1
            if pas > 0
                J(1:pas,:)	= j(length(J)-pas+1:length(J),:);
                J(pas+1:length(J),:)	= j(1:length(J)-pas,:);
            elseif pas < 0
                J(1:length(J)+pas,:)	= j(-pas+1:length(J),:);
                J(length(J)+pas+1:length(J),:)	= j(1:-pas,:);
            end
        elseif dim == 2
            if pas > 0
                J(:,1:pas)	= j(:,length(J)-pas+1:length(J));
                J(:,pas+1:length(J))	= j(:,1:length(J)-pas);
            elseif pas < 0
                J(:,1:length(J)+pas)	= j(:,-pas+1:length(J));
                J(:,length(J)+pas+1:length(J))	= j(:,1:-pas);
            end
        end
    elseif ndims(J) == 3
        if dim == 1
            if pas > 0
                J(1:pas,:,:)	= j(length(J)-pas+1:length(J),:,:);
                J(pas+1:length(J),:,:)	= j(1:length(J)-pas,:,:);
            elseif pas < 0
                J(1:length(J)+pas,:,:)	= j(-pas+1:length(J),:,:);
                J(length(J)+pas+1:length(J),:,:)	= j(1:-pas,:,:);
            end
        elseif dim == 2
            if pas > 0
                J(:,1:pas,:)	= j(:,length(J)-pas+1:length(J),:);
                J(:,pas+1:length(J),:)	= j(:,1:length(J)-pas,:);
            elseif pas < 0
                J(:,1:length(J)+pas,:)	= j(:,-pas+1:length(J),:);
                J(:,length(J)+pas+1:length(J),:)	= j(:,1:-pas,:);
            end
        elseif dim == 3
            if pas > 0
                J(:,:,1:pas)	= j(:,:,length(J)-pas+1:length(J));
                J(:,:,pas+1:length(J))	= j(:,:,1:length(J)-pas);
            elseif pas < 0
                J(:,:,1:length(J)+pas)	= j(:,:,-pas+1:length(J));
                J(:,:,length(J)+pas+1:length(J))	= j(:,:,1:-pas);
            end
        end
    end
end

