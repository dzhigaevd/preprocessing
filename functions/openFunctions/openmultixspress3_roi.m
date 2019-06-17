function out = openmultixspress3_roi(varargin)
% --- USAGE
% --- OUT = OPENMULTIEIGER4M(FILENAME)        -- direct call to 
% --- open all images of a HDF5 file
% --- OUT = OPENMULTIEIGER4M(FILENAME,STARTINDEX,ENDINDEX)-- call to 
% --- open a block of images (lines). 
% --- OUT = OPENMULTIEIGER4M(FILENAME,STARTINDEX,ENDINDEX,ROI)-- call to 
% --- open an ROI of a block of images. 
% --- OUT = OPENMULTIEIGER4M(FILENAME,ROI)-- call to 
% --- open an ROI of all images. 
% ---
% --- Input Argument
% --- FILE        image file name
% --- STARTINDEX  starting image index
% --- ENDINDEX    ending image index
% --- ROI         [roiYstart, roiYend, roiXstart, roiXend]
% --- ROIextended [roiYstart, roiYend, roiXstart, roiXend, frameStart, frameEnd]

% --- Output Argument
% --- OUT  structure containing header & image block.
% ---
% --- D. Dzhigaev
% --- $Revision 1.0 $Date 20190411 $ function to open merlin files
% --- Adaptation from P10, PETRA III scripts by M. Sprung

out = [];       

% =========================================================================
% --- distinguish case used and create variable file
% =========================================================================
if ( nargin == 1 )
    file               = varargin{1}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiStart          = 1                                                  ;        
    roiRange          = inf                                                ;
elseif ( nargin == 2 )
    file               = varargin{1}                                       ;
    StartIndex         = 1                                                 ;
    EndIndex           = inf                                               ;
    roiStart           = varargin{2}(1)                                    ;
    roiEnd             = varargin{2}(2)                                    ;    
    roiRange           = roiEnd - roiStart + 1                             ;
    clear roiEnd                                                           ;
elseif ( nargin == 3 )
    file               = varargin{1}                                       ;
    StartIndex         = varargin{2}                                       ;
    EndIndex           = varargin{3}                                       ;
    if StartIndex >= EndIndex
        dummy      = EndIndex                                              ;
        EndIndex   = StartIndex                                            ;
        StartIndex = dummy                                                 ;
        clear dummy                                                        ;
    end
    roiStart          = 1                                                 ;    
    roiRange          = inf                                               ;    
elseif ( nargin == 4 )
    file               = varargin{1}                                       ;
    StartIndex         = varargin{2}                                       ;
    EndIndex           = varargin{3}                                       ;
    if StartIndex >= EndIndex
        dummy      = EndIndex                                              ;
        EndIndex   = StartIndex                                            ;
        StartIndex = dummy                                                 ;
        clear dummy                                                        ;
    end
    roiStart          = varargin{4}(1)                                    ;
    roiEnd            = varargin{4}(2)                                    ;        
    frameStart        = varargin{4}(3)                                    ;
    frameEnd          = varargin{4}(4)                                    ;
    
    roiRange          = roiEnd - roiStart + 1                              ;
    
    clear roiEnd                                                           ;
end
    
    
% =========================================================================
% --- open file
% =========================================================================
[fid,message] = fopen(file,'r')                                            ;        
if ( fid == -1 )                                                             % return if open fails
    uiwait(msgbox(message,'File Open Error','error','modal'))              ;
    return                                                                 ;
end
fullhfile          = fopen(fid)                                            ;  % get full pathname
[pathstr,name,ext] = fileparts(fullhfile)                                  ;
fclose(fid)                                                                ;

% =========================================================================
% --- check if 'file' points to a master or a data file
% =========================================================================
master = 0                                                                 ; % switch for master file
data   = 0                                                                 ; % switch for data file
data1  = 0                                                                 ; % switch for first data file
if strcmp(ext,'.h5')
    master     = 1                                                         ;
    masterfile = file                                                      ;
    disp('Using master file!');
elseif strcmp(ext,'.hdf5')
    data = 1                                                               ;
    datafile = file;
    disp('Using data file!');
else
    disp('Error! Input File not recognized as master or data file')        ;
    return
end

if data == 1
    fileInfo = h5info(file)                                                ;
    dataInfo = h5info(file, '/entry_0000/measurement/Merlin/data/')        ;
    
    % dataFullSize: [xDiff,yDiff,fastAxis,slowAxis]
    dataFullSize  = [dataInfo.Dataspace.Size, length(fileInfo.Groups)]     ;
    
    if EndIndex == inf
        EndIndex = dataInfo.Dataspace.Size(3)                              ; % set EndIndex to the chunk size
    end
    
    if dataInfo.Dataspace.Size(3) >= StartIndex && dataInfo.Dataspace.Size(3) >= EndIndex                    
        try
            out = zeros([dataFullSize(1),dataFullSize(2),frameEnd - frameStart,EndIndex - StartIndex],'single') ;
        catch EM
            disp(EM)                                               ;
            fprintf('Please load the data in smaller chunks! \n')     ;
            return                                                 ;
        end
        
        for ii = 1 : EndIndex - StartIndex + 1
            for jj = 1 : frameEnd - frameStart + 1                
                out(:,:,jj,ii) = (h5read(file,sprintf('/entry_%04i/measurement/Merlin/data/',StartIndex+ii-1),[roiXstart roiStart frameStart+jj-1],[roiXrange roiRange 1]))                                 ;                
            end
            fprintf('Read line #%d \n',ii);
        end
        
        fprintf('Final data size [%d,%d,%d,%d] \n', (size(out)));
    else
        fprintf('Error! Data file contains only %i frames!',dataInfo.Dataspace.Size(3)) ; %DSPS
        return
    end 
end


end

