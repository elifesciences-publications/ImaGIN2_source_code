function[data,h]=ImaGIN_SONGetADCChannel(fid,chan,ReadTime)
% Reads an ADC (waveform) channel from a SON file.
%
%

% Malcolm Lidierth 02/02

% Modified 15/12/03
% Now returns header correctly if frames are of
% unequal length.  Also, FrameLength calculation changed to reflect unequal
% lengths
%
% Bug Fix 15/12/04
% h.stop now set to end of last block when NumFrames=1 (previously set to
% end of first block).
%
% Bug Fix 28/05/04
% Revision of 15/12/04 fix. h.stop now set correctly when NumFrames==1 and
% there are fewer than 5 data blocks.

Info=SONChannelInfo(fid,chan);
if(Info.kind==0) 
    warning('ImaGIN_SONGetADCChannel: No data on that channel');
    return;
end;

FileH=SONFileHeader(fid);
SizeOfHeader=20;                                            % Block header is 20 bytes long
header=SONGetBlockHeaders(fid,chan);

NumberOfSamples=sum(header(5,:));                           % Sum of samples in all blocks
SampleInterval=(header(3,1)-header(2,1))/(header(5,1)-1);   % Sample interval in clock ticks

%OD: find blocks to read
TimeSamples=cumsum(header(5,:));
ReadStart=max(find(TimeSamples<=ReadTime.Start));
if isempty(ReadStart)
    ReadStart=1;
end
ReadStop=min(find(TimeSamples>=ReadTime.Stop));
if isempty(ReadStop)
    ReadStop=Info.blocks;
end
ReadStop=min([Info.blocks ReadStop]);
NumberOfSamples=sum(header(5,ReadStart:ReadStop));

if(nargout>1)
h.FileName=Info.FileName;                                   % Set up the header information to return
h.system=['SON' num2str(FileH.systemID)];                   % if wanted
h.FileChannel=chan;
h.phyChan=Info.phyChan;
h.kind=Info.kind;
%h.blocks=Info.blocks;
%h.preTrig=Info.preTrig;
h.comment=Info.comment;
h.title=Info.title;
h.sampleinterval=SONGetSampleInterval(fid,chan);
h.scale=Info.scale;
h.offset=Info.offset;
h.units=Info.units;
end;


NumFrames=1;                                                % Number of frames. Initialize to one.
Frame(1:ReadStart)=1;
% for i=1:Info.blocks-1       %OD                                % Check for discontinuities in data record
for i=ReadStart:ReadStop-1                                    
    IntervalBetweenBlocks=header(2,i+1)-header(3,i);
    if IntervalBetweenBlocks>SampleInterval                 % If true data is discontinuous (triggered)
        NumFrames=NumFrames+1;                              % Count discontinuities (NumFrames)
        Frame(i+1)=NumFrames;                               % Record the frame number that each block belongs to
    else
        Frame(i+1)=Frame(i);                                % Pad between discontinuities
    end;
end;

if NumFrames==1                                             % Continuous sampling - one frame only
    data=int16(zeros(1,NumberOfSamples));                   % Pre-allocate memory for data
    pointer=1;
    h.start=header(2,1);                                    % Time of first data point (clock ticks)
    [nrows,ncols]=size(header);
    h.stop=header(3,ncols);                                 % End of data (clock ticks) Bug fix 28/5/05
    
%     for i=1:Info.blocks       %OD                                  
    for i=ReadStart:ReadStop                                        
        fseek(fid,header(1,i)+SizeOfHeader,'bof');
        data(pointer:pointer+header(5,i)-1)=fread(fid,header(5,i),'int16=>int16');
        pointer=pointer+header(5,i);
    end;
    %OD
    h.start=header(2,ReadStart);                                    % Time of first data point (clock ticks)
    h.stop=header(3,ReadStop);                                 % End of data (clock ticks)

else                                                        % Frame based data -  multiple frames
    FrameLength=max(histc(Frame,[1:NumFrames]))*max(header(5,:));% Maximum data points to a frame
    data=int16(zeros(NumFrames,FrameLength));               % Pre-allocate array    
    start=1;                                                % Pointer into array for each disk data block
    Frame(Info.blocks+1)=-99;                               % Dummy entry to avoid index error in for loop
    
    h.start(1)=header(2,1);                  % Time of first data point in frame #1 (clock ticks)****************15/12/03
%     for i=1:Info.blocks       %OD                                  
    for i=ReadStart:ReadStop                                        
        fseek(fid,header(1,i)+SizeOfHeader,'bof');
        data(Frame(i),start:start+header(5,i)-1)=fread(fid,header(5,i),'int16=>int16');
        if Frame(i+1)==Frame(i)
            start=start+header(5,i);                        % Increment pointer or.....
        else
            start=1;                                        % begin new frame
            h.stop(Frame(i))=header(3,i);                   % End time for this frame, clock ticks
            if(i<Info.blocks)                               %***************15/12/03
                h.start(Frame(i+1))=header(2,i+1);              % Time of first data point in next frame (clock ticks)*************15/12/03
            end;
        end;
    end;
end;

h.start=SONTicksToSeconds(fid,h.start);
h.stop=SONTicksToSeconds(fid,h.stop);
