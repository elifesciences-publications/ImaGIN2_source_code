function J	= ImaGIN_shift(J,pas,dim)

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

